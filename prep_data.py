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
COUNTRY_VECTOR_PATH = 'data/countries_iso3_md5_6fb2431e911401992e6e56ddf0a9bcda.gpkg'
EEZ_VECTOR_PATH = 'data/eez_iso_sov1_md5_d18f061b8628dc6da36067db7b485d3a.gpkg'
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
        LOGGER.debug(f'{vector_path} {field_id}')
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


def main():
    """Entry point."""
    task_graph = taskgraph.TaskGraph('.', multiprocessing.cpu_count(), 15.0)
    worker_list = []
    target_path_set = set()
    stitch_raster_task_list = []
    solution_lookup = collections.defaultdict(list)
    for solution_dir in glob.glob(os.path.join(SOLUTIONS_DIR, '*')):
        if not os.path.isdir(solution_dir):
            continue
        raster_step_set = collections.defaultdict(list)
        raster_list = list(glob.glob(os.path.join(solution_dir, '*.tif')))
        step_set = set()
        for raster_path in raster_list:
            match = re.match('.*[^\d](\d+)\.tif', raster_path)
            if match:
                percent_fill = int(match.group(1))
                solution_lookup[os.path.basename(solution_dir)].append(
                    percent_fill)
                raster_step_set[percent_fill].append(raster_path)
        for percent_fill, raster_path_list in raster_step_set.items():
            scenario_id = os.path.basename(solution_dir)
            target_raster_path = os.path.join(
                STICH_DIR, f'{scenario_id}_{percent_fill}.tif')
            merge_task = task_graph.add_task(
                func=merge_raster_set,
                args=(raster_path_list, target_raster_path),
                target_path_list=[target_raster_path],
                task_name=f'merge for {target_raster_path}')
            stitch_raster_task_list.append((target_raster_path, scenario_id, percent_fill, merge_task))

    eez_ids_names = task_graph.add_task(
        func=get_field_names,
        args=(EEZ_VECTOR_PATH, EEZ_FIELD_ID), store_result=True).get()
    country_ids_names = task_graph.add_task(
        func=get_field_names,
        args=(COUNTRY_VECTOR_PATH, COUNTRY_FIELD_ID), store_result=True).get()

    rasterized_dict = collections.defaultdict(dict)
    scenario_percent_to_hash = dict()
    for stitch_raster_path, scenario_id, percent_fill, dependent_task in stitch_raster_task_list:
        stitch_raster_info = task_graph.add_task(
            func=pygeoprocessing.get_raster_info,
            args=(stitch_raster_path,),
            dependent_task_list=[dependent_task],
            store_result=True,
            task_name=f'get raster info for {stitch_raster_path}')
        stitch_hash = hashlib.sha1(pickle.dumps(stitch_raster_info.get())).hexdigest()
        scenario_percent_to_hash[(scenario_id, percent_fill)] = stitch_hash
        if stitch_hash not in rasterized_dict:
            for prefix, vector_path, vector_field_id, field_names in [
                    ('eez', EEZ_VECTOR_PATH, EEZ_FIELD_ID, eez_ids_names),
                    ('country', COUNTRY_VECTOR_PATH, COUNTRY_FIELD_ID,
                     country_ids_names),
                    ]:
                global_target_path = os.path.join(
                    MASK_DIR, f'{prefix}_global_{stitch_hash}.tif')
                rasterize_task = task_graph.add_task(
                    func=rasterize_with_base,
                    args=(vector_path, stitch_hash, global_target_path),
                    target_path_list=[global_target_path],
                    task_name=f'rasterize {global_target_path}')
                rasterized_dict[stitch_hash][f'{prefix}_global'] = (
                    global_target_path, rasterize_task)
                for field_name in field_names:
                    local_target_path = os.path.join(
                        MASK_DIR, f'{prefix}_{field_name}_{stitch_hash}.tif')
                    rasterize_task = task_graph.add_task(
                        func=rasterize_with_base,
                        args=(vector_path, stitch_hash, local_target_path),
                        kwargs={
                            'field_id': vector_field_id,
                            'field_value': field_name},
                        target_path_list=[local_target_path],
                        task_name=f'rasterize {local_target_path}')
                    rasterized_dict[stitch_hash][f'{prefix}_{field_name}'] = (
                        local_target_path, rasterize_task)

    # for each scenario
    scenario_to_stats_map = collections.defaultdict(dict)
    for (optimization_mask_raster_path, scenario_id, percent_fill, merge_task) in \
            stitch_raster_task_list:
        # for global & each country
        stitch_hash = scenario_percent_to_hash[(scenario_id, percent_fill)]
        scenario_stats_task = task_graph.add_task(
            func=calculate_pixel_stats,
            args=(
                optimization_mask_raster_path, scenario_id, percent_fill,
                stitch_hash, country_ids_names, rasterized_dict),
            dependent_task_list=[merge_task],
            store_result=True,
            task_name=f'stats for {scenario_id}{percent_fill}')
        scenario_to_stats_map[(scenario_id, percent_fill)] = \
            scenario_stats_task

    # stitch_hash = hashlib.sha1(pickle.dumps(stitch_raster_info))
    # to get a task do this: rasterized_dict[stitch_hash][f'{prefix}_{field_name}']
    task_graph.close()
    task_graph.join()
    task_graph = None
    LOGGER.info(f'all done: {scenario_to_stats_map}')


def calculate_pixel_stats(
        optimization_mask_raster_path, scenario_id, percent_fill, stitch_hash,
        country_ids_names, country_eez_mask_dict):
    """Calc all pixel counts in country/eez total and those chosen in opt.

    Args:
        optimization_mask_raster_path (str): path to optimization scenario mask
        scenario_id (str): overall scenario id like "A" or "L1".
        percent_fill (int/str): sub scenario showing what % was filled up.
        stitch_hash (str): used to index into country_eez_mask_dict to get the
            correctly sized masks
        country_ids_names (list): list of country id names.
        country_eez_mask_dict (str): indexed by stitch and then
            '{eez/country}_{country/eez code}' and returns a mask that is
            uniquely that.

    Returns:
        {
            'n_pixels_land': x,
            'n_pixels_eez': x,
            'n_pixels_total': x,
            'sqkm_land': x,
            'sqkm_eez': x,
            'sqkm_total': x,
            'n_pixels_land_selected': x,
            'n_pixels_eez_selected': x,
            'n_pixels_total_selected': x,
            'sqkm_land_selected': x,
            'sqkm_eez_selected': x,
            'sqkm_total_selected': x,
        }
    """
    # TODO: this
    for country_id in ['global'] + country_ids_names:
        mask_raster_path, mask_rasterize_task = \
            country_eez_mask_dict[stitch_hash][f'country_{country_id}']
        LOGGER.debug(
            f'optimize stats for: {optimization_mask_raster_path} '
            f'{scenario_id}:{percent_fill}:stitch_hash:{country_id}')


if __name__ == '__main__':
    main()
