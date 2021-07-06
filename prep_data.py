"""Entrypoint script to process and stitch data."""
import collections
import glob
import logging
import os
import re
import multiprocessing
import tempfile
import shutil

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
logging.getLogger('taskgraph').setLevel(logging.INFO)
LOGGER = logging.getLogger(__name__)


SOLUTIONS_DIR = 'data/solutions'
STICH_DIR = 'churn/target_stitch_dir'
COUNTRY_VECTOR_PATH = 'data/countries_iso3_md5_6fb2431e911401992e6e56ddf0a9bcda.gpkg'
EEZ_VECTOR_PATH = 'data/eez_iso_sov1_md5_d18f061b8628dc6da36067db7b485d3a.gpkg'
EEZ_FIELD_ID = 'ISO_SOV1'
COUNTRY_FIELD_ID = 'iso3'
RES_KM = 2.0
os.makedirs(STICH_DIR, exist_ok=True)


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
    n_cols, n_rows = base_raster_info['raster_size']

    raster_driver = gdal.GetDriverByName(
        DEFAULT_GTIFF_CREATION_TUPLE_OPTIONS[0])
    new_raster = raster_driver.Create(
        target_raster_path, n_cols, n_rows, 1, base_raster_info['datatype'],
        options=DEFAULT_GTIFF_CREATION_TUPLE_OPTIONS[1])
    new_raster.SetProjection(base_raster_info['projection_wkt'])
    new_raster.SetGeoTransform(base_raster_info['geotransform'])
    new_band = new_raster.GetRasterBand(1)
    new_band.SetNoDataValue(base_raster_info['nodata'][0])
    new_band = None
    new_raster = None

    where_clause = None
    if field_id is not None:
        where_clause = f'"{field_id}"="{field_value}"'
    pygeoprocessing.rasterize(
        vector_path, target_raster_path, burn_values=(1,), option_list=None,
        where_clause=where_clause)


def _zonal_stats(base_raster_path, vector_path):
    """Calculate zonal stats of base over vector."""
    raster_info = pygeoprocessing.get_raster_info(base_raster_path)
    base_dir = os.path.join(os.path.dirname(vector_path), 'tmp')
    os.makedirs(base_dir, exist_ok=True)
    working_dir = tempfile.mkdtemp(dir=base_dir)
    reprojected_vector_path = os.path.join(
        working_dir, os.path.basename(vector_path))
    pygeoprocessing.reproject_vector(
        vector_path, raster_info['projection_wkt'],
        reprojected_vector_path,
        driver_name='GPKG', copy_fields=True)
    stats = pygeoprocessing.zonal_statistics(
        (base_raster_path, 1), reprojected_vector_path,
        polygons_might_overlap=False)
    shutil.rmtree(working_dir)
    return stats


def main():
    """Entry point."""
    task_graph = taskgraph.TaskGraph('.', multiprocessing.cpu_count()//2, 15.0)
    stitch_raster_task_list = []
    scenario_percent_type_map = collections.defaultdict(
        lambda: collections.defaultdict(dict))
    for solution_dir in glob.glob(os.path.join(SOLUTIONS_DIR, '*')):
        if not os.path.isdir(solution_dir):
            continue

        raster_step_set = collections.defaultdict(list)
        raster_list = list(glob.glob(os.path.join(solution_dir, '*.tif')))
        for raster_path in raster_list:
            match = re.match('.*[^\d](\d+)\.tif', raster_path)
            percent_fill = int(match.group(1))
            raster_step_set[percent_fill].append(raster_path)
        for percent_fill, raster_path_list in raster_step_set.items():
            scenario_id = os.path.basename(solution_dir)
            merged_raster_path = os.path.join(
                STICH_DIR, f'{scenario_id}_{percent_fill}.tif')
            merge_task = task_graph.add_task(
                func=merge_raster_set,
                args=(raster_path_list, merged_raster_path),
                target_path_list=[merged_raster_path],
                task_name=f'merge for {merged_raster_path}')
            stitch_raster_task_list.append(
                (merged_raster_path, scenario_id, percent_fill, merge_task))
            LOGGER.debug(f'{merged_raster_path}, {EEZ_VECTOR_PATH}')
            eez_stats_task = task_graph.add_task(
                func=_zonal_stats,
                args=(merged_raster_path, EEZ_VECTOR_PATH),
                dependent_task_list=[merge_task],
                store_result=True,
                task_name=f'eez stats for {merged_raster_path}')
            country_stats_task = task_graph.add_task(
                func=_zonal_stats,
                args=(merged_raster_path, COUNTRY_VECTOR_PATH),
                dependent_task_list=[merge_task],
                store_result=True,
                task_name=f'country stats for {merged_raster_path}')

            scenario_percent_type_map[scenario_id][percent_fill]['eez'] = \
                eez_stats_task
            scenario_percent_type_map[scenario_id][percent_fill]['country'] = \
                country_stats_task

    vector_fid_field_map = collections.defaultdict(
        lambda: collections.defaultdict(lambda: None))
    for vector_path, vector_id, vector_field in [
            (EEZ_VECTOR_PATH, 'eez', EEZ_FIELD_ID),
            (COUNTRY_VECTOR_PATH, 'country', COUNTRY_FIELD_ID)]:
        vector = gdal.OpenEx(vector_path, gdal.OF_VECTOR)
        layer = vector.GetLayer()
        for feature in layer:
            vector_fid_field_map[vector_id][feature.GetField(vector_field)] = \
                feature.GetFID()

    table_file = open('result.csv', 'w')
    table_file.write(
        'scenario,country,res_km,n_pu,n_selected,n_pu_land,n_selected_land,'
        'prop_selected,sqkm_selected,prop_selected_land,sqkm_selected_land\n')
    for scenario_id, percent_region_type_map in \
            scenario_percent_type_map.items():
        for percent_fill, region_to_stats_map in \
                percent_region_type_map.items():
            eez_stats = region_to_stats_map["eez"].get()
            country_stats = region_to_stats_map["country"].get()
            table_file.write(
                f'{scenario_id},global')
            for country_id, country_fid in \
                    vector_fid_field_map['country'].items():
                eez_fid = vector_fid_field_map['eez'][country_id]
                LOGGER.debug(
                    f'{scenario_id},{percent_fill},{country_id},{country_fid}:\n'
                    f'\t{country_stats[country_fid]}\n')
                if eez_fid is not None:
                    LOGGER.debug(f'\t{eez_stats[eez_fid]}\n')
                    pass
                table_file.write(
                    f'''{scenario_id},{country_id},{get_stats(
                        country_fid, RES_KM, country_stats[country_fid],
                        eez_stats[country_fid])}\n''')
    table_file.close()
    task_graph.close()
    task_graph.join()
    task_graph = None
    LOGGER.info('all done')


def get_stats(country_fid, res_km, country_stats_dict, eez_stats_dict):
    """Return formatted stats."""
    n_selected_land = country_stats_dict[country_fid]['count']
    n_pu_land = country_stats_dict[country_fid]['nodata_count'] + n_selected_land
    prop_selected_land = n_selected_land / n_pu_land

    n_selected = n_selected_land
    n_pu = n_pu_land
    if country_fid in eez_stats_dict:
        n_selected_eez = eez_stats_dict[country_fid]['count']
        n_pu_eez = eez_stats_dict[country_fid]['nodata_count'] + n_selected_eez
        n_selected += n_selected_eez
        n_pu += n_pu_eez
    prop_selected = n_selected / n_pu
    sqkm_selected = n_selected * res_km
    sqkm_selected_land = n_selected_land * res_km
    return f'{res_km},{n_pu},{n_selected},{n_pu_land},{n_selected_land},{prop_selected},{sqkm_selected},{prop_selected_land},{sqkm_selected_land}'


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
