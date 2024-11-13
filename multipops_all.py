import pandas as pd
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.mask import mask
from shapely.geometry import shape
from rasterio.features import geometry_mask
from rasterstats import zonal_stats
import matplotlib.pyplot as plt


import pdb # pdb.set_trace()

import sys
import os
from damIA_impact import population_impact, increase_resolution, download_url


### Input parameters
reservoirs_path='file3.geojson'
minarea = 1 # km2
targetres = 10  # [m]
startyear = 1975
endyear = 2020
downloaddir = '/download/'

### Preparations
gdf_in = gpd.read_file(reservoirs_path)

refyears = np.arange(startyear, endyear+1, 5)
bins = list(refyears[1:])  # 5-9 years offset from Year of Completion
labels = list(refyears[:-2])
#bins = list(refyears[2:])   # 10-14 years offset from Year of Completion
#labels = list(refyears[:-3])
labelsnp = np.asarray(labels)

# Start year filter
gdf_past_styear = gdf_in.loc[(gdf_in['Year of Completion'] >= bins[0])].copy()

# Minimum area filter (according to polygon size, since ICOLD area data is extremely noisy)
reservoirs = gdf_past_styear.loc[(gdf_past_styear['plg_a_km2'] >= 1.0)].copy()

# Exclude transboundary dams to avoid "diluting" information
transb_dams = ['YACYRETA','TAIPINGWAN', 'ITAIPU', 'MANGLA'] #(TAIPINGWAN is transboundary China/N.Korea, ITAIPU Paraguay/Brazil, MANGLA Pakistan/India (also in Worldpop))
reservoirs = reservoirs.loc[~(reservoirs['Name of the dam'].isin(transb_dams))]

# Bin dams into the map reference years (i.e. 5-9 years before Year of Completion)
reservoirs['refyear'] = pd.cut(x = reservoirs['Year of Completion'], bins=bins, labels=labels, include_lowest=True, right=False).astype(int)

pdb.set_trace()

## Input raster file paths
popgrid_names = ['GHS-POP 2023','LandScan','WorldPop','WorldPop UN-adjusted','GRUMP','GRUMP UN-adjusted','GWP','GWP UN-adjusted']
popgrid_shortnames = ['ghs','lscan','wpop','wpop-un','grump','grump-un','gwp','gwp-un']
#popgrid_names = ['GRUMPv1','GWPv4.11','LandScan']
#popgrid_shortnames = ['grump','gwp','lscan']
srcres = [100, 1000, 100, 100, 1000, 1000, 1000, 1000] # Native resolutions of popgrids [m]. Manually provided since some popgrids are in wgs84 and others not
basepath = '/'
# for each popdens source, make a list of paths to all years in refyears (None if no map exists)
ghs_files = ['api_ghs','api_ghs','api_ghs','api_ghs','api_ghs','api_ghs','api_ghs','api_ghs']
#ghs1k_files = [basepath+'/GHSL2023/GHS_POP_E1975_GLOBE_R2023A_4326_30ss_V1_0/GHS_POP_E1975_GLOBE_R2023A_4326_30ss_V1_0.tif',basepath+'/GHSL2023/GHS_POP_E1980_GLOBE_R2023A_4326_30ss_V1_0/GHS_POP_E1980_GLOBE_R2023A_4326_30ss_V1_0.tif',basepath+'/GHSL2023/GHS_POP_E1985_GLOBE_R2023A_4326_30ss_V1_0/GHS_POP_E1985_GLOBE_R2023A_4326_30ss_V1_0.tif', basepath+'/GHSL2023/GHS_POP_E1990_GLOBE_R2023A_4326_30ss_V1_0/GHS_POP_E1990_GLOBE_R2023A_4326_30ss_V1_0.tif', basepath+'/GHSL2023/GHS_POP_E1995_GLOBE_R2023A_4326_30ss_V1_0/GHS_POP_E1995_GLOBE_R2023A_4326_30ss_V1_0.tif', basepath+'/GHSL2023/GHS_POP_E2000_GLOBE_R2023A_4326_30ss_V1_0/GHS_POP_E2000_GLOBE_R2023A_4326_30ss_V1_0.tif', basepath+'/GHSL2023/GHS_POP_E2005_GLOBE_R2023A_4326_30ss_V1_0/GHS_POP_E2005_GLOBE_R2023A_4326_30ss_V1_0.tif', basepath+'/GHSL2023/GHS_POP_E2010_GLOBE_R2023A_4326_30ss_V1_0/GHS_POP_E2010_GLOBE_R2023A_4326_30ss_V1_0.tif']
wp_files = [None,None,None,None,None,'api_wp','api_wp','api_wp']
wpun_files = [None,None,None,None,None,'api_wpun','api_wpun','api_wpun']
grump_files = [None,None,None,basepath+'grumpv1/gl_grumpv1_pcount_90_ascii_30/glup90g.tif',basepath+'grumpv1/gl_grumpv1_pcount_95_ascii_30/glup95g.tif',basepath+'grumpv1/gl_grumpv1_pcount_00_ascii_30/glup00g.tif',None,None]
grumpun_files = [None,None,None,basepath+'grumpv1/gl_grumpv1_pcount_90_ascii_30/glup90ag.tif',basepath+'grumpv1/gl_grumpv1_pcount_95_ascii_30/glup95ag.tif',basepath+'grumpv1/gl_grumpv1_pcount_00_ascii_30/glup00ag.tif',None,None]
gwp_files = [None,None,None,None,None,basepath+'/GWPv411/gpw-v4-population-count-rev11_2000_30_sec_tif/gpw_v4_population_count_rev11_2000_30_sec.tif',basepath+'GWPv411/gpw-v4-population-count-rev11_2005_30_sec_tif/gpw_v4_population_count_rev11_2005_30_sec.tif',basepath+'GWPv411/gpw-v4-population-count-rev11_2010_30_sec_tif/gpw_v4_population_count_rev11_2010_30_sec.tif']
gwpun_files = [None,None,None,None,None,basepath+'/GWPv411/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11_2000_30_sec_tif/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2000_30_sec.tif',basepath+'GWPv411/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11_2005_30_sec_tif/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2005_30_sec.tif',basepath+'GWPv411/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11_2010_30_sec_tif/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2010_30_sec.tif']
lscan_files = [None,None,None,None,None,basepath+'Landscan/landscan-global-2000-assets/landscan-global-2000.tif',basepath+'Landscan/landscan-global-2005-assets/landscan-global-2005.tif',basepath+'Landscan/landscan-global-2010-assets/landscan-global-2010.tif']
raster_files = np.transpose([ghs_files, lscan_files, wp_files, wpun_files, grump_files, grumpun_files, gwp_files, gwpun_files])

#pdb.set_trace()

### Main program

def main(reservoirs, raster_files, popgrid_shortnames):

    # Iterate over each polygon
    for idx, reservoir in reservoirs.iterrows():

        print("")
        #print("*** Processing reservoir "+str(idx+1)+"/"+str(len(reservoirs))+': '+reservoir['Name of the dam'])
        print("*** Processing reservoir: " +str(idx) + ' ' + reservoir['Name of the dam'])

        # Task 0: List of paths of popdens grids corresponding to dam reference year
        map_paths = raster_files[np.where(labelsnp == reservoir.refyear)]

        # Loop over different population data sources for the reference year
        for id_src, map in enumerate(map_paths[0]):
            if map is not None:
                if map == 'api_ghs':    # Use damIA to download and compute GHS-POP
                    pop_sum, population_map_year = population_impact(reservoirs, idx, reservoir.refyear, downloaddir) # Compute ghs population using damIA
                else:
                    if map == 'api_wp':     # Download country data from WorldPop (does not work for transboundary reservoirs)
                        iso3 = reservoir.ISO_CODES
                        wp_url = 'https://data.worldpop.org/GIS/Population/Global_2000_2020/' + str(reservoir.refyear) + '/' + iso3 + '/' + iso3.lower() + '_ppp_' + str(reservoir.refyear) + '.tif'
                        download_filepath = basepath + 'WorldPop/' + str(reservoir.refyear) + '_unconstr/' + iso3.lower() + '_ppp_' + str(reservoir.refyear) + '.tif'
                        map = worldpop_download(wp_url, download_filepath)
                    if map == 'api_wpun':     # Download UN-adjusted country data from WorldPop (does not work for transboundary reservoirs)
                        iso3 = reservoir.ISO_CODES
                        wpun_url = 'https://data.worldpop.org/GIS/Population/Global_2000_2020/' + str(reservoir.refyear) + '/' + iso3 + '/' + iso3.lower() + '_ppp_' + str(reservoir.refyear) + '_UNadj.tif'
                        download_filepath = basepath + 'WorldPop/' + str(reservoir.refyear) + '_unconstr_UNadj/' + iso3.lower() + '_ppp_' + str(reservoir.refyear) + '_UNadj.tif'
                        map = worldpop_download(wpun_url, download_filepath)

                    #Crop raster to the extent of the polygon
                    cropped_raster, cropped_meta = crop_raster_to_polygon(map, reservoir)

                    #Increase resolution of the cropped raster
                    hres_raster, hres_meta = increase_resolution(cropped_raster, cropped_meta, srcres[id_src], targetres=targetres)
                    #with rasterio.open('output.tif', 'w', **hres_meta) as dst:     # testing
                    #    dst.write(hres_raster[0,:,:], 1)

                    # Set nodata cells to zero
                    hres_raster[np.where(hres_raster<0)] = 0

                    # Calculate the pixel sum for high-resolution pixels inside the polygon
                    stats = zonal_stats(reservoir.geometry, hres_raster, affine=hres_meta["transform"], stats="sum")
                    pop_sum = stats[0]["sum"]

                # Add the pixel sum as a new attribute to the polygon
                reservoirs.at[idx, popgrid_shortnames[id_src]] = pop_sum

    pdb.set_trace()

    # Save the updated polygons GeoJSON file
    #output_file = "multipops_ICOLD2023.geojson"
    output_file = "multipops_ICOLD2023_1014y.geojson"
    reservoirs.to_file(output_file, driver="GeoJSON")

    pdb.set_trace()

### Helper functions

def crop_raster_to_polygon(raster_file, polygon):
    with rasterio.open(raster_file) as src:
        out_image, out_transform = mask(src, [polygon.geometry], all_touched=True, crop=True, pad=True)
        out_meta = src.meta.copy()

        if len(out_image.shape) > 2:
            out_image = out_image[0,:,:]
        # with rasterio.open('output.tif', 'w', **hres_meta) as dst:     # testing
        #    dst.write(hres_raster[0,:,:], 1)

    out_meta.update({"driver": "GTiff", "height": out_image.shape[0], "width": out_image.shape[1], "transform": out_transform})
    return out_image, out_meta


def worldpop_download(wp_url, download_filepath):
    if not os.path.exists(download_filepath):
        print("        *** DOWNLOADING WORLDPOP POPULATION DENSITY MAP ***")
        print('        URL: ', wp_url)
        download_url(wp_url, download_path=download_filepath)
    map = download_filepath
    return map

# def increase_resolution(raster, in_meta, srcres, targetres=10):
#     factor = srcres / targetres
#     new_shape = (raster.shape[0], raster.shape[1] * factor, raster.shape[2] * factor)
#     hres_raster = np.zeros(new_shape)
#
#     for band in range(raster.shape[0]):
#         for i in range(raster.shape[1]):
#             for j in range(raster.shape[2]):
#                 hres_raster[band, i * factor : (i + 1) * factor, j * factor : (j + 1) * factor] = raster[band, i, j] / (factor**2)    # population count is divided into the high-resolution cells
#
#     out_meta = in_meta.copy()
#     out_transform = in_meta['transform']
#     out_meta.update({"driver": "GTiff", "height": hres_raster.shape[1], "width": hres_raster.shape[2], "transform": out_transform})
#     out_meta.update(transform=out_transform * out_transform.scale(1/factor, 1/factor))
#
#     return hres_raster, out_meta

### Execute main program
main(reservoirs, raster_files, popgrid_shortnames)
