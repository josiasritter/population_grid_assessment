
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
# MAIN PROGRAM FOR ESTIMATING DISPLACEMENTS

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
