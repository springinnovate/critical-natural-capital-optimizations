"""Entrypoint script to process and stitch data."""
import hashlib
import pickle
import collections
import glob
import logging
import os
import re
import multiprocessing

from pygeoprocessing.geoprocessing_core import DEFAULT_GTIFF_CREATION_TUPLE_OPTIONS
from osgeo import gdal
import pygeoprocessing
import taskgraph

gdal.SetCacheMax(2**27)

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'))
logging.getLogger('taskgraph').setLevel(logging.WARN)
LOGGER = logging.getLogger(__name__)


SOLUTIONS_DIR = 'data/solutions'
STICH_DIR = 'churn/target_stitch_dir'
EEZ_VECTOR_PATH = 'data/countries_iso3_md5_6fb2431e911401992e6e56ddf0a9bcda.gpkg'
COUNTRY_VECTOR_PATH = 'data/eez_iso_sov1_md5_d18f061b8628dc6da36067db7b485d3a.gpkg'
MASK_DIR = 'churn/mask_dir'
EEZ_FIELD_ID = 'ISO_SOV1'
COUNTRY_FIELD_ID = 'iso3'

os.makedirs(STICH_DIR, exist_ok=True)
os.makedirs(MASK_DIR, exist_ok=True)


def main():
    """Entry point."""
    pass


def merge_raster_set(raster_path_list, target_raster_path):
    """Merge all the rasters into one target.

    Args:
        raster_path_list (list): a list of raster paths.

    Return:
        None
    """
    projection_set = set()
    pixel_size_set = set()
    bounding_box_list = []
    for raster_path in raster_path_list:
        raster_info = pygeoprocessing.get_raster_info(raster_path)
        projection_set.add(raster_info['projection_wkt'])
        bounding_box_list.append(raster_info['bounding_box'])
        pixel_size_set.add(raster_info['pixel_size'])
    target_bounding_box = pygeoprocessing.merge_bounding_box_list(
        bounding_box_list, 'union')
    pixel_size = next(iter(pixel_size_set))

    n_cols = int(abs(
        (target_bounding_box[2]-target_bounding_box[0])/pixel_size[0]))+1
    n_rows = int(abs(
        (target_bounding_box[1]-target_bounding_box[3])/pixel_size[1]))+1
    print(pixel_size_set, n_cols, n_rows, projection_set, target_bounding_box)

    raster_driver = gdal.GetDriverByName(
        DEFAULT_GTIFF_CREATION_TUPLE_OPTIONS[0])
    new_raster = raster_driver.Create(
        target_raster_path, n_cols, n_rows, 1, raster_info['datatype'],
        options=DEFAULT_GTIFF_CREATION_TUPLE_OPTIONS[1])
    new_raster.SetProjection(next(iter(projection_set)))
    new_raster.SetGeoTransform(
        [target_bounding_box[0], pixel_size[0], 0.0,
         target_bounding_box[3], 0.0, pixel_size[1]])
    new_band = new_raster.GetRasterBand(1)
    new_band.SetNoDataValue(raster_info['nodata'][0])
    new_band = None
    new_raster = None
    pygeoprocessing.stitch_rasters(
        [(path, 1) for path in raster_path_list],
        ['near']*len(raster_path_list),
        (target_raster_path, 1),
        overlap_algorithm='etch',
        area_weight_m2_to_wgs84=False)


def get_field_names(vector_path, field_id):
    """Return set of all field names found with vector_path and field_id."""
    vector = gdal.OpenEx(vector_path, gdal.OF_VECTOR)
    layer = vector.GetLayer()
    field_name_set = set()
    for feature in layer:
        field_name_set.add(feature.GetField(field_id))
    return field_name_set


def rasterize_with_base(
        vector_path, base_raster_info, target_raster_path, field_id=None,
        field_value=None):
    """Rasterize vector onto target with the same size as base.

    Args:
        vector_path (str): path to vector to rasterize
        base_raster_info (dict): dictionary to make base raster target
            size.
        target_raster_path (str): path to desired target
        field_id/field_value (str/str): if not none, restrict the polygon
            rasterization to these fields.

    Returns:
        None
    """
    pass


if __name__ == '__main__':
    task_graph = taskgraph.TaskGraph('.', multiprocessing.cpu_count(), 15.0)
    worker_list = []
    target_path_set = set()
    stitch_raster_task_list = []
    for solution_dir in glob.glob(os.path.join(SOLUTIONS_DIR, '*')):
        if not os.path.isdir(solution_dir):
            continue
        raster_step_set = collections.defaultdict(list)
        raster_list = list(glob.glob(os.path.join(solution_dir, '*.tif')))
        step_set = set()
        for raster_path in raster_list:
            match = re.match('.*[^\d](\d+)\.tif', raster_path)
            if match:
                raster_step_set[int(match.group(1))].append(raster_path)
        for step_size, raster_path_list in raster_step_set.items():
            print(step_size)
            target_raster_path = os.path.join(
                STICH_DIR, f'{os.path.basename(solution_dir)}_{step_size}.tif')
            merge_task = task_graph.add_task(
                func=merge_raster_set,
                args=(raster_path_list, target_raster_path),
                target_path_list=[target_raster_path])
            stitch_raster_task_list.append((target_raster_path, merge_task))

    eez_ids_names = task_graph.add_task(
        func=get_field_names,
        args=(EEZ_VECTOR_PATH, EEZ_FIELD_ID)).get()
    country_ids_names = task_graph.add_task(
        func=get_field_names,
        args=(COUNTRY_VECTOR_PATH, COUNTRY_FIELD_ID)).get()
    #TODO: don't forget to check that these are the same

    rasterized_dict = collections.defaultdict(dict)
    for stitch_raster_path, dependent_task in stitch_raster_task_list:
        stitch_raster_info = task_graph.add_task(
            func=pygeoprocessing.get_raster_info,
            args=(stitch_raster_path,),
            dependent_task_list=[dependent_task])
        stitch_hash = hashlib.sha1(pickle.dumps(stitch_raster_info))
        if stitch_hash not in rasterized_dict:
            for prefix, vector_path, vector_field_id, field_names in [
                    ('eez', EEZ_VECTOR_PATH, EEZ_FIELD_ID, eez_ids_names),
                    ('country', COUNTRY_VECTOR_PATH, COUNTRY_FIELD_ID,
                     country_ids_names),
                    ]:
                global_target_path = os.path.join(
                    MASK_DIR, f'{prefix}_global_{stitch_hash}.tif')
                task_graph.add_task(
                    func=rasterize_with_base,
                    args=(vector_path, stitch_hash, target_raster_path),
                    target_path_list=[target_raster_path],
                    task_name=f'rasterize {target_raster_path}')
                rasterized_dict[stitch_hash][f'{prefix}_global'] = \
                    target_raster_path
                for field_name in field_names:
                    local_target_path = os.path.join(
                        MASK_DIR, f'{prefix}_{field_name}_{stitch_hash}.tif')
                    rasterized_dict[stitch_hash][f'{prefix}_{field_name}'] = \
                        target_raster_path
                    task_graph.add_task(
                        func=rasterize_with_base,
                        args=(vector_path, stitch_hash, target_raster_path),
                        kwargs={
                            'field_id': vector_field_id,
                            'field_value': field_name},
                        target_path_list=[target_raster_path],
                        task_name=f'rasterize {target_raster_path}')

    # stitch_hash = hashlib.sha1(pickle.dumps(stitch_raster_info))
    # to get a task do this: rasterized_dict[stitch_hash][f'{prefix}_{field_name}']
    task_graph.close()
    task_graph.join()
    task_graph = None
