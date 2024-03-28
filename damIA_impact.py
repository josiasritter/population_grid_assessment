
# IMPORT LIBRARIES -----------------------------------------------------------
import math
import os
import numpy as np

import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterstats import zonal_stats
import shapely.speedups

import time
from datetime import date
import calendar
import requests
from itertools import product
from zipfile import ZipFile

from damIA_main import create_directories, cell_area


# (REMOVE AFTER TESTING)
# import matplotlib.pyplot as plt
import pdb # pdb.set_trace()

#######################################################################################################################
# MAIN PROGRAM FOR ESTIMATING IMPACTS


"""
#TESTING

from damIA_impact import *
run_name = 'Stimson2000'
resshape_path = '/Users/ritterj1/PythonProjects/PycharmProjects/damIA/reservoir_performance/stimson_validationdams.geojson'
impact(resshape_path, run_name=run_name, dam_name='unnamed_dam', commission_year=2000, population='on', landcover='on', protection='on', indigenous='on')
# REMEMBER THAT I OVERWRITE commission_year FOR THIS RUN IN LINE 180! UNDO THIS FOR OTHER RUNS!

from damIA_impact import *
run_name = 'CGIAR_present2'
resshape_path = '/Users/ritterj1/PythonProjects/PycharmProjects/damIA/CGIAR_present2/alldams_CGIAR_present2.geojson'
impact(resshape_path, run_name=run_name, dam_name='unnamed_dam', commission_year=2000, population='on', landcover='on', protection='on', indigenous='on')
# REMEMBER THAT I OVERWRITE commission_year FOR THIS RUN IN LINE 180! UNDO THIS FOR OTHER RUNS!

from damIA_impact import *
run_name = 'GeoDAR_ICOLD_displaced_valid'
resshape_path = '/Users/ritterj1/PythonProjects/PycharmProjects/damIA/GeoDAR_ICOLD_displaced_valid/output/GeoDAR_p1975_displ.geojson'
impact(resshape_path, run_name=run_name, dam_name='unnamed_dam', commission_year=2023, population='on', landcover='off', protection='off', indigenous='off')

from damIA_impact import *
run_name = 'GeoDAR_test1'
resshape_path = '/Users/ritterj1/PythonProjects/PycharmProjects/damIA/GeoDAR_test3/output/geodar_test3.geojson'
impact(resshape_path, run_name=run_name, dam_name='unnamed_dam', commission_year=2015, population='on', landcover='on', protection='on', indigenous='on')

resshape_path = '/Users/ritterj1/PythonProjects/dam_IA/damIA_v4/Tuoba/output/Tuoba.geojson' 
resshape = fiona.open(resshape_path, 'r')
respoly = shape(resshape[0]['geometry'])
downloaddir='/Users/ritterj1/PythonProjects/dam_IA/damIA_v4/Tuoba/download/'
outdir='/Users/ritterj1/PythonProjects/dam_IA/damIA_v4/Tuoba/output/'
commission_year=2022
lost_protected_sites, lost_protected_area = protection_impact(respoly, downloaddir, outdir)
# landcover_impact(resshape, commission_year, downloaddir)

"""

def impact(resshape_path, run_name='unnamed_run', dam_name='unnamed_dam',
           commission_year=2023, population='off', landcover='off', protection='off', indigenous='off', silentmode='on'):

    if silentmode == 'on':  # Temporarily suppress all Python warnings in command line
        import warnings
        warnings.filterwarnings('ignore')

    # Check input keywords
    if population == 'off' and landcover == 'off' and protection == 'off' and indigenous == 'off':
        raise Exception("Input Error: All impact categories are switched off. At least one of the keywords < population >, < landcover >, or < protection > must have the value 'on' to compute impacts.")

    # Create download and output directories
    downloaddir, outdir, catchdir, dam_name_clean = create_directories(run_name, dam_name)

    # Load reservoir polygon file and metadata
    reservoir_gdf = gpd.read_file(resshape_path)
    n_reservoirs = len(reservoir_gdf)

    # Initialise lists of impacts to be filled in the loop below
    displaced_population_allres, population_map_year_allres, lost_built_area_allres, lost_crop_area_allres, lost_forest_area_allres, landcover_map_year_allres, lost_protected_sites_allres, lost_protected_area_allres, lost_indigenous_area_allres = [], [], [], [], [], [], [], [], []

    # Loop over reservoirs in polygon file
    for i in range(n_reservoirs):
        print("")
        print("COMPUTING IMPACTS FOR DAM " + str(i+1) + ' of ' + str(n_reservoirs))
        #pdb.set_trace()

        # Excecution time. Remove after testing!
        st = time.time()

        # Check for commission year
        #commission_year = reservoir_gdf['Year of Completion'][i] # GeoDAR inputs %%%
        #commission_year = reservoir_gdf['COD'][i]  # CGIAR inputs %%%
        if commission_year < 1975 and population == 'on':
            raise Exception("No exposure data available for such early years! This tool provides reliable estimates of displaced people only for dams completed after around 1980.")
        elif commission_year < 1980 and population == 'on':
            print("Warning! Earliest population map is from 1975. This tool provides reliable estimates of displaced population only for dams completed after around 1980.")
        # elif commission_year < 2005 and landcover == 'on':
            # print('Warning! Land cover map is from year 2000. This tool provides reliable estimates of displaced land covers only for dams completed after around 2005.')

        # Displaced population
        if population == 'on':
            print("    *** DISPLACED PEOPLE ***")
            #if i == 0:
            #   reservoir_gdf_esri54009 = reservoir_gdf.to_crs('ESRI:54009')    # Change crs of reservoir polygons to crs of population density data
            #displaced_population, population_map_year = population_impact(reservoir_gdf_esri54009, i, commission_year, downloaddir)
            displaced_population, population_map_year = population_impact(reservoir_gdf, i, commission_year, downloaddir)
            displaced_population_allres.append(displaced_population)
            population_map_year_allres.append(population_map_year)

            if i == n_reservoirs-1:     # Write estimated impacts to output geodataframe
                reservoir_gdf['displaced_pop'] = displaced_population_allres
                reservoir_gdf['pop_map_year'] = population_map_year_allres

        # Lost land covers
        if landcover == 'on':
            print("    *** LOST LAND COVERS ***")
            lost_built_area, lost_crop_area, lost_forest_area, landcover_map_year = landcover_impact(reservoir_gdf, i, commission_year, downloaddir)
            lost_built_area_allres.append(lost_built_area)
            lost_crop_area_allres.append(lost_crop_area)
            lost_forest_area_allres.append(lost_forest_area)
            landcover_map_year_allres.append(landcover_map_year)

            if i == n_reservoirs-1:     # Write estimated impacts to output geodataframe
                reservoir_gdf['lost_built_km2'] = lost_built_area_allres
                reservoir_gdf['lost_crop_km2'] = lost_crop_area_allres
                reservoir_gdf['lost_forest_km2'] = lost_forest_area_allres
                reservoir_gdf['landc_map_year'] = landcover_map_year_allres

        # Lost environmentally and culturally protected sites and areas. Works only for future dams! No past analysis, as no past data available.
        if protection == 'on':
            print("    *** LOST PROTECTED SITES AND AREAS ***")
            if i == 0:
                # Load data of protected sites and areas into memory
                sites_gdf, areas_gdf = prepare_protection_data(downloaddir)
                shapely.speedups.enable()

            respoly = reservoir_gdf.loc[i, 'geometry']
            lost_protected_sites, lost_protected_area = protection_impact(respoly, sites_gdf, areas_gdf, outdir)
            lost_protected_sites_allres.append(lost_protected_sites)
            lost_protected_area_allres.append(lost_protected_area)

            if i == n_reservoirs-1:     # Write estimated impacts to output geodataframe
                reservoir_gdf['lost_prot_sites'] = lost_protected_sites_allres
                reservoir_gdf['lost_prot_km2'] = lost_protected_area_allres

        # Lost indigenous lands. Only for future dams! No past analysis, as no past data available.
        if indigenous == 'on':
            print("    *** LOST INDIGENOUS PEOPLES' LANDS ***")
            if i == 0:
                # Load indigenous peoples' lands data into memory (data is not available online but upon individual request to first author of https://www.nature.com/articles/s41893-018-0100-6)
                indigenous_shp = '/Users/ritterj1/GIS_data/IPL_IndigenousPeoplesLands_2017/01_Data/IPL_IndigenousPeoplesLands_2017/IPL_2017_wgs84.shp'
                indigenous_gdf = gpd.read_file(indigenous_shp)
                shapely.speedups.enable()

            respoly = reservoir_gdf.loc[i, 'geometry']
            lost_indigenous_area = indigenous_impact(respoly, indigenous_gdf, outdir)
            lost_indigenous_area_allres.append(lost_indigenous_area)

            if i == n_reservoirs-1:     # Write estimated impacts to output geodataframe
                reservoir_gdf['lost_indig_km2'] = lost_indigenous_area_allres

        # Excecution times. Remove after testing!
        et = time.time()
        print('Execution time: ', str(et - st), ' seconds')

    # Overwrite reservoir polygon file to save impact estimates
    reservoir_gdf.to_file(resshape_path)

    #pdb.set_trace()


#######################################################################################################################
# PROGRAM COMPONENTS (IMPACT CATEGORIES)

def population_impact(reservoir_gdf, i, commission_year, downloaddir):

    # Select reference year of population map (maps available: for 1975-2030 in 5y intervals; https://ghsl.jrc.ec.europa.eu/download.php?ds=pop)
    #reference_year = max(1975, 5 * math.floor((commission_year - 8.6) / 5))   # "Large dams on average take 8.6 years to build (std = 4.6 y)" (https://www.sciencedirect.com/science/article/pii/S0301421513010926); round reference year down to nearest 5 year interval of population density maps.
    reference_year = commission_year  # Fit to classification in popgrid study
    #reference_year = 2000 # CGIAR present2 run %%%
    reference_year = str(reference_year)


    # Select population data tile(s)
    bounds = reservoir_gdf['geometry'][i].bounds    # bounds of reservoir
    tilenames = create_tilenames_pop(bounds)
    #url_base = 'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2022A/GHS_POP_E'
    url_base = 'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E'

    # Create download directory
    #population_dir = downloaddir + 'population' + reference_year + '/'
    population_dir = downloaddir + 'population_ghsl2023_' + reference_year + '/'
    if not os.path.exists(population_dir):
        os.mkdir(population_dir)

    # Set count of population to zero to prepare for adding population from several data tiles
    displaced_population_alltiles = 0

    for tile in tilenames:
        # Download population density data
        # url = url_base + reference_year + '_GLOBE_R2022A_54009_100/V1-0/tiles/GHS_POP_E' + reference_year + '_GLOBE_R2022A_54009_100_V1_0_' + tile + '.zip'
        # url = url_base + reference_year + '_GLOBE_R2023A_54009_100/V1-0/tiles/GHS_POP_E' + reference_year + '_GLOBE_R2023A_54009_100_V1_0_' + tile + '.zip'  # Mollweide
        url = url_base + reference_year + '_GLOBE_R2023A_4326_3ss/V1-0/tiles/GHS_POP_E' + reference_year + '_GLOBE_R2023A_4326_3ss_V1_0_' + tile + '.zip' # WGS84
        tile_outpath = population_dir + tile + '.zip'
        #tile_tifpath = population_dir + 'GHS_POP_E' + reference_year + '_GLOBE_R2022A_54009_100_V1_0_' + tile + '.tif'
        #tile_tifpath = population_dir + 'GHS_POP_E' + reference_year + '_GLOBE_R2023A_54009_100_V1_0_' + tile + '.tif'
        tile_tifpath = population_dir + 'GHS_POP_E' + reference_year + '_GLOBE_R2023A_4326_3ss_V1_0_' + tile + '.tif'

        if not os.path.exists(tile_tifpath):
            print("        *** DOWNLOADING POPULATION DENSITY MAP TILE ***")
            print('        URL: ', url)
            download_url(url, download_path=tile_outpath)

            # Unzip and remove zipfile
            with ZipFile(tile_outpath) as zObject:
                zObject.extractall(path=population_dir)
            os.remove(tile_outpath)
            os.remove(population_dir + 'GHSL_Data_Package_2023_light.pdf')
            os.remove(population_dir + 'GHS_POP_GLOBE_R2023A_input_metadata.xlsx')

        # Crop population density map to reservoir area
        with rasterio.open(tile_tifpath) as population_grid:
            population_cropped, population_cropped_transform = mask(population_grid, [reservoir_gdf['geometry'][i]], crop=True, pad=True, nodata=0)  # maybe needs another .mask
            population_cropped = population_cropped[0]
            population_cropped_meta = population_grid.meta.copy()
            population_cropped_meta.update({"driver": "GTiff", "height": population_cropped.shape[0], "width": population_cropped.shape[1], "transform": population_cropped_transform})

        # testing
        # outras = rasterio.open('geodar_test1_pol1_displ' + '.tif', 'w', driver='GTiff', height=population_cropped.shape[0], width=population_cropped.shape[1], count=1, dtype=population_cropped.dtype, nodata=0, crs='+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ', transform=population_cropped_transform)
        # outras.write(population_cropped, 1)
        # outras.close()

        # Increase resolution to 10 m before applying zonal stats (to reduce uncertainties related to location of cell centres)
        hres_raster, out_meta = increase_resolution(population_cropped, population_cropped_meta, population_cropped_transform[0], targetres=population_cropped_transform[0]/10)

        # Set nodata cells to zero
        hres_raster[np.where(hres_raster < 0)] = 0

        displaced_population = zonal_stats(reservoir_gdf['geometry'][i], hres_raster, affine=out_meta['transform'], stats=['sum'])
        if bool(displaced_population[0]['sum']):
            #displaced_population_sum = np.around(displaced_population[0]['sum'], decimals=0)
            displaced_population_sum = np.around(displaced_population[0]['sum'], decimals=3)
        else:
            displaced_population_sum = 0
        displaced_population_alltiles += displaced_population_sum

    # Print results
    print("            *** Overall number of displaced people: ", displaced_population_alltiles)

    return displaced_population_alltiles, int(reference_year)


def landcover_impact(reservoir_gdf, i, commission_year, downloaddir):

    # Select reference year of land cover map (two maps available: year 2000 and year 2020; https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/download.html)
    if commission_year <= 2023:
        reference_year = '2000' # land cover map of year 2000 selected, because in map of 2020 land cover is likely already modified by dam construction
    else: reference_year = '2020'
    if commission_year < 2000:  # If comission year is earlier than land cover map, exit the land cover impact module and return NaN as outputs
        reference_year, lost_built_alltiles, lost_crop_alltiles, lost_forest_alltiles = np.nan, np.nan, np.nan, np.nan  # float("nan"), float("nan"), float("nan"), float("nan")
        print("            *** Warning: No lost land cover computed because dam commission year is before the earliest available land cover map (year 2000)")
        return reference_year, lost_built_alltiles, lost_crop_alltiles, lost_forest_alltiles

    # Select land cover data tile(s)
    bounds = reservoir_gdf['geometry'][i].bounds  # bounds of reservoir
    tilenames = create_tilenames_lc(bounds)
    url_base = 'https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/'

    # Create download directory
    landcover_dir = downloaddir + 'landcover' + reference_year + '/'
    if not os.path.exists(landcover_dir):
        os.mkdir(landcover_dir)

    # Set counts of lost land cover by type to zero to prepare for adding lost land cover from several data tiles
    lost_built_alltiles, lost_crop_alltiles, lost_forest_alltiles = 0, 0, 0

    for tile in tilenames:
        # Download land cover data
        url = url_base + reference_year + '/' + tile + '.tif'
        tile_outpath = landcover_dir + tile + '.tif'
        if not os.path.exists(tile_outpath):
            print("        *** DOWNLOADING LAND COVER MAP TILE ***")
            print('        URL: ', url)
            download_url(url, download_path=tile_outpath)

        # Crop land cover to reservoir area
        with rasterio.open(tile_outpath) as landcover_grid:
            landcover_cropped, landcover_cropped_transform = mask(landcover_grid, [reservoir_gdf['geometry'][i]], crop=True)
            landcover_cropped = landcover_cropped[0]

        # Create cell area grid for land cover in reservoir area
        areagrid = cell_area(landcover_cropped, landcover_cropped_transform)

        # Create one grid per land cover class with the cell area as values
        built_codes = [250]
        crop_codes = [244]
        forest_codes = list(range(25,49)) + list(range(125,149))    # first range: forest on terra firme with tree height > 3m; second range: forest in wetlands with tree height >3m

        lost_built_tile = compute_lost_landcover_by_type(built_codes, landcover_cropped, areagrid, [reservoir_gdf['geometry'][i]], landcover_cropped_transform)
        lost_crop_tile = compute_lost_landcover_by_type(crop_codes, landcover_cropped, areagrid, [reservoir_gdf['geometry'][i]], landcover_cropped_transform)
        lost_forest_tile = compute_lost_landcover_by_type(forest_codes, landcover_cropped, areagrid, [reservoir_gdf['geometry'][i]], landcover_cropped_transform)

        lost_built_alltiles += lost_built_tile
        lost_crop_alltiles += lost_crop_tile
        lost_forest_alltiles += lost_forest_tile

    # Print results
    print("            *** Overall size of lost built area: ", lost_built_alltiles, "km^2")
    print("            *** Overall size of lost cropland area: ", lost_crop_alltiles, "km^2")
    print("            *** Overall size of lost forest area: ", lost_forest_alltiles, "km^2")

    return lost_built_alltiles, lost_crop_alltiles, lost_forest_alltiles, int(reference_year)


def protection_impact(respoly, sites_gdf, areas_gdf, outdir):

    # Protected sites (points): find those inside the reservoir polygon
    sites_gdf['geometry'] = sites_gdf.representative_point()  # Convert multipoint geometries to point geometries
    boolean_mask = sites_gdf.within(respoly)  # Filter for points inside reservoir polygon # THIS NEEDS TESTING!
    lost_sites_gdf = sites_gdf.loc[boolean_mask]

    if len(lost_sites_gdf) > 0:     # Count lost protected sites and write their details to hard disk
        lost_protected_sites = len(lost_sites_gdf)
        lost_sites_gdf.to_file(outdir + 'lost_protected_sites.geojson')     # THIS NEEDS TO BE ADJUSTED TO BATCH PROCESSING! CURRENTLY, OUTPUT FILES ARE CREATED FOR EACH DAM AND PROGRESSIVELY OVERWRITTEN!
    else: lost_protected_sites = 0

    # Protected areas (polygons): Find and clip those intersecting with the reservoir polygon
    polys = areas_gdf.geometry
    intersects = polys.intersection(respoly)
    boolean_mask = ~intersects.is_empty
    lost_areas_gdf = areas_gdf.loc[boolean_mask].copy()
    lost_areas_gdf['geometry'] = intersects[boolean_mask]

    if len(lost_areas_gdf) > 0:     # Calculate overall lost protected area and write the areas and their details to hard disk
        lost_areas_gdf_mollweide = lost_areas_gdf.to_crs('ESRI:54009').copy()  # Reproject to Mollweide projection (since area calculation does not work with unprojected coordinates)
        lost_areas_gdf['lost_area_km2'] = lost_areas_gdf_mollweide['geometry'].area / 1e6
        lost_protected_area = np.around(lost_areas_gdf['lost_area_km2'].sum(), decimals=4)
        lost_areas_gdf.to_file(outdir + 'lost_protected_areas.geojson')     # THIS NEEDS TO BE ADJUSTED TO BATCH PROCESSING! CURRENTLY, OUTPUT FILES ARE CREATED FOR EACH DAM AND PROGRESSIVELY OVERWRITTEN!
    else: lost_protected_area = 0

    # Print results
    print("            *** Overall number of lost protected sites: ", lost_protected_sites)
    if lost_protected_sites > 0:
        print("            *** For details, see output file: lost_protected_sites.geojson")
    print("            *** Overall size of lost protected area: ", lost_protected_area, "km^2")
    if lost_protected_area > 0:
        print("            *** For details, see output file: lost_protected_areas.geojson")

    return lost_protected_sites, lost_protected_area


def indigenous_impact(respoly, indigenous_gdf, outdir):

    # Find and clip indigenous poeples' areas intersecting with reservoir polygon
    polys = indigenous_gdf.geometry
    intersects = polys.intersection(respoly)
    boolean_mask = ~intersects.is_empty
    lost_areas_gdf = indigenous_gdf.loc[boolean_mask].copy()
    lost_areas_gdf['geometry'] = intersects[boolean_mask]

    # Write lost indigenous area polygons and their details to hard disk, and calculate overall lost indigenous area
    if len(lost_areas_gdf) > 0:
        lost_areas_mollweide = lost_areas_gdf.to_crs('ESRI:54009').copy()  # Reproject to Mollweide projection for area calculation (does not work with unprojected coordinates)
        lost_areas_gdf['lost_area_km2'] = lost_areas_mollweide['geometry'].area / 1000000
        lost_areas_gdf.to_file(outdir + 'lost_indigenous_areas.geojson')        # THIS NEEDS TO BE ADJUSTED TO BATCH PROCESSING! CURRENTLY, OUTPUT FILES ARE CREATED FOR EACH DAM AND PROGRESSIVELY OVERWRITTEN!
        lost_indigenous_area = np.around(lost_areas_gdf['lost_area_km2'].sum(), decimals=4)
    else: lost_indigenous_area = 0

    # Print results
    print("            *** Overall size of lost indigenous area: ", lost_indigenous_area, "km^2")
    if lost_indigenous_area > 0:
        print("            *** For details, see output file: lost_indigenous_areas.geojson")

    return lost_indigenous_area


#######################################################################################################################
# HELPER FUNCTIONS

# WGS84
def create_tilenames_pop(bounds):
    """For given bounds of a reservoir polygon, create the names of the required population density data tiles (see https://ghsl.jrc.ec.europa.eu/download.php?ds=pop)"""
    round_bounds = [math.floor((bounds[0] + 180.008) / 10) + 1, math.floor((89.100 - bounds[1]) / 10) + 1, math.floor((bounds[2] + 180.008) / 10) + 1, math.floor((89.100 - bounds[3]) / 10) + 1]  # Round bounds to population density data tile top left corner.

    colnrs = np.unique([round_bounds[0], round_bounds[2]])
    colnames = ['C' + str(colnr) for colnr in colnrs]
    rownrs = np.unique([round_bounds[1], round_bounds[3]])
    rownames = ['R' + str(rownr) for rownr in rownrs]

    # Create list of combinations of tile indices in case more than one tile is needed
    xy_combinations = list(list(zip(rownames, element)) for element in product(colnames, repeat=len(rownames)))
    xy_combs_flat = [y for x in xy_combinations for y in x]

    # Names of download tiles
    xy_combs_str = [combination[0] + '_' + combination[1] for combination in xy_combs_flat]
    tilenames = list(np.unique(xy_combs_str))

    return tilenames

# # Mollweide
# def create_tilenames_pop(bounds):
#     """For given bounds of a reservoir polygon, create the names of the required population density data tiles (see https://ghsl.jrc.ec.europa.eu/download.php?ds=pop)"""
#     round_bounds = [math.floor((bounds[0] + 18.041e6) / 1e6) + 1, math.floor((9e6 - bounds[1]) / 1e6) + 1, math.floor((bounds[2] + 18.041e6) / 1e6) + 1, math.floor((9e6 - bounds[3]) / 1e6) + 1]  # Round bounds to population density data tile top left corner.
#
#     colnrs = np.unique([round_bounds[0], round_bounds[2]])
#     colnames = ['C' + str(colnr) for colnr in colnrs]
#     rownrs = np.unique([round_bounds[1], round_bounds[3]])
#     rownames = ['R' + str(rownr) for rownr in rownrs]
#
#     # Create list of combinations of tile indices in case more than one tile is needed
#     xy_combinations = list(list(zip(rownames, element)) for element in product(colnames, repeat=len(rownames)))
#     xy_combs_flat = [y for x in xy_combinations for y in x]
#
#     # Names of download tiles
#     xy_combs_str = [combination[0] + '_' + combination[1] for combination in xy_combs_flat]
#     tilenames = list(np.unique(xy_combs_str))
#
#     return tilenames


def create_tilenames_lc(bounds):
    """For given bounds of a reservoir polygon, create the names of the required land cover data tiles (see https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/download.html)"""
    round_bounds = [math.floor(bounds[0]/10)*10, math.ceil(bounds[1]/10)*10, math.floor(bounds[2]/10)*10, math.ceil(bounds[3]/10)*10]    # Round bounds to land cover data tile top left corner (10 degrees x 10 degrees tiles).
    lons = np.unique([round_bounds[0], round_bounds[2]])
    lats = np.unique([round_bounds[1], round_bounds[3]])

    # Convert lons and lats from negative/positive values to directional values
    lons_WE = []
    for lon in lons:
        if lon <= -100:
            lons_WE.append(str(abs(lon)) + 'W')
        elif lon < 0:
            lons_WE.append('0'+ str(abs(lon)) + 'W')
        elif lon < 100:
            lons_WE.append('0'+ str(lon) + 'E')
        elif lon >= 100:
            lons_WE.append(str(lon) + 'E')
    lats_SN = []
    for lat in lats:
        if lat < 0:
            lats_SN.append(str(abs(lat)) + 'S')
        else:
            lats_SN.append(str(lat) + 'N')

    # Create list of combinations of tile indices in case more than one tile is needed
    lonlat_combinations = list(list(zip(lats_SN, element)) for element in product(lons_WE, repeat=len(lats_SN)))
    lonlat_combs_flat = [y for x in lonlat_combinations for y in x]

    # Names of download tiles
    lonlat_combs_str = [combination[0] + '_' + combination[1] for combination in lonlat_combs_flat]
    tilenames = list(np.unique(lonlat_combs_str))

    return tilenames


def download_url(url, download_path="./"):
    """Download data from a given url to a given download path on the hard drive"""
    while True:
        try:
            response = requests.get(url, timeout=30)
        except:
            print("    *** WARNING! Download froze, likely due to network instability. Restarting download...")
            time.sleep(2)
        else:
            if response.ok:
                # inmemory = tiff.imread(io.BytesIO(response.content))  # reads download data directly to memory (without writing to disk)
                with open(download_path, 'wb') as savefile:
                    savefile.write(response.content)  # write download data to disk
                break

    # while True:
    #     try:
    #         response = requests.get(url, timeout=10)
    #         if response.ok:
    #             # inmemory = tiff.imread(io.BytesIO(response.content))  # reads download data directly to memory (without writing to disk)
    #             with open(download_path, 'wb') as savefile:
    #                 savefile.write(response.content)  # write download data to disk
    #             break
    #         else:
    #             Exception("Download link failed: ", url)
    #     except requests.exceptions.Timeout as err:
    #         print(err)
    #     time.sleep(2)
    #


def increase_resolution(raster, in_meta, srcres, targetres=10):
    factor = round(srcres / targetres)
    new_shape = (raster.shape[0] * factor, raster.shape[1] * factor)
    hres_raster = np.zeros(new_shape)

    #for band in range(raster.shape[0]):
    for i in range(raster.shape[0]):
        for j in range(raster.shape[1]):
            hres_raster[i * factor: (i + 1) * factor, j * factor: (j + 1) * factor] = raster[i, j] / (factor ** 2)  # population count is divided into the high-resolution cells

    out_meta = in_meta.copy()
    out_transform = in_meta['transform']
    out_meta.update({"driver": "GTiff", "height": hres_raster.shape[0], "width": hres_raster.shape[1], "transform": out_transform})
    out_meta.update(transform=out_transform * out_transform.scale(1/factor, 1/factor))

    return hres_raster, out_meta


def prepare_protection_data(downloaddir):
    # Create download directory
    protection_dir = downloaddir + 'protection/'
    if not os.path.exists(protection_dir):
        os.mkdir(protection_dir)

    # Download protected areas if not yet present
    if len(os.listdir(protection_dir)) == 0:
        download_protection_data(protection_dir)

    # Create list of shapefiles, as the downloaded dataset consists of 3 point files (protected sites) and 3 polygon files (protected areas)
    pointfiles, polyfiles = [], []
    for root, d_names, f_names in os.walk(protection_dir):
        for f in f_names:
            if f.endswith('points.shp'):
                pointfiles.append(os.path.join(root, f))
            if f.endswith('polygons.shp'):
                polyfiles.append(os.path.join(root, f))

    if not bool(pointfiles) or not bool(polyfiles):
        raise Exception("WARNING! Data of protected sites and areas is incomplete!!! Remove folder <protection> in download directory and restart program!")

    for i, file in enumerate(pointfiles):
        if i == 0:
            sites_gdf = gpd.read_file(file)
        else:
            add_sites_gdf = gpd.read_file(file)
            sites_gdf = pd.concat([sites_gdf, add_sites_gdf])

    for i, file in enumerate(polyfiles):
        if i == 0:
            areas_gdf = gpd.read_file(file)
        else:
            add_areas_gdf = gpd.read_file(file)
            areas_gdf = pd.concat([areas_gdf, add_areas_gdf])

    return sites_gdf, areas_gdf


def download_protection_data(protection_dir):
    # Prepare download url and path
    url_base = 'https://d1gam3xoknrgr2.cloudfront.net/current/'  # URL is automatically updated according to current month, since data is updated and renamed monthly: https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA
    today_monthyear = calendar.month_abbr[date.today().month] + str(date.today().year)
    url_filename = 'WDPA_' + today_monthyear + '_Public_shp.zip'
    url = url_base + url_filename
    download_path = protection_dir + '/' + url_filename

    # Start download
    print("        *** DOWNLOADING PROTECTED AREAS AND SITES. This may take several minutes (file > 1 GB) ***")
    print('        URL: ', url)
    download_url(url, download_path=download_path)
    print("        *** UNZIPPING PROTECTED AREAS AND SITES ***")
    unzip_all(dir=protection_dir)


def unzip_all(dir='.'):
    """Iteratively unzip files in given directory and subdirectories until no zip files remain"""
    zipfiles = True
    while bool(zipfiles):
        zipfiles = []
        for root, d_names, f_names in os.walk(dir):
            for f in f_names:
                if f.endswith('.zip'):
                    zipfiles.append(os.path.join(root, f))
        if bool(zipfiles):
            for file in zipfiles:
                #targetdir = dir + '/' + file[:-4]
                targetdir = file[:-4]
                if not os.path.exists(targetdir):
                    os.mkdir(targetdir)
                with ZipFile(file, 'r') as zObject:
                    zObject.extractall(path=targetdir)
                os.remove(file)


def compute_lost_landcover_by_type(landcover_codes, landcover_grid, areagrid, respoly, affine):
    """Sum the area with given list of landcover codes inside the reservoir polygon"""
    lost_areagrid = np.zeros(np.shape(landcover_grid))
    landcover_cells = np.where(np.isin(landcover_grid, landcover_codes))
    lost_areagrid[landcover_cells] = areagrid[landcover_cells]
    #outras = rasterio.open('geodar_test1_pol1_lostcrop' + '.tif', 'w', driver='GTiff', height=lost_areagrid.shape[0], width=lost_areagrid.shape[1], count=1, dtype=lost_areagrid.dtype, nodata=-9999, crs='+proj=latlong', transform=affine)
    #outras.write(lost_areagrid, 1)
    lost_landcover = zonal_stats(respoly, lost_areagrid, affine=affine, stats=['sum'])
    if bool(lost_landcover[0]['sum']):
        lost_landcover_sum = np.around(lost_landcover[0]['sum'], decimals=4)
    else:
        lost_landcover_sum = 0

    return lost_landcover_sum
