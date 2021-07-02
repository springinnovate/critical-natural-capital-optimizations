"""Entrypoint script to process and stitch data."""
import collections
import glob
import logging
import os
import re

from pygeoprocessing.geoprocessing_core import DEFAULT_GTIFF_CREATION_TUPLE_OPTIONS
from osgeo import gdal
import pygeoprocessing

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'))
logging.getLogger('taskgraph').setLevel(logging.WARN)
LOGGER = logging.getLogger(__name__)


SOLUTIONS_DIR = 'data/solutions'
STICH_DIR = 'target_stitch_dir'
os.makedirs(STICH_DIR, exist_ok=True)


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
    print(target_raster_path)
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


if __name__ == '__main__':
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
            merge_raster_set(raster_path_list, target_raster_path)
            break
        break
