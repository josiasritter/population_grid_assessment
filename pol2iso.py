import geopandas as gpd
import pdb # pdb.set_trace()

def get_country_iso_codes(geojson_file, countries_shapefile):
    data = gpd.read_file(geojson_file)
    countries = gpd.read_file(countries_shapefile)
    pdb.set_trace()

    iso_codes_list = []
    for geometry in data['geometry']:
        countries_intersected = countries[countries.intersects(geometry)]
        iso_codes = countries_intersected['ISO3_CODE'].tolist()
        iso_codes_str = ','.join(iso_codes)
        iso_codes_list.append(iso_codes_str)

    data['ISO_CODES'] = iso_codes_list  # Add ISO_CODES column to the GeoDataFrame

    return data

# Execute
#geojson_file='/Users/ritterj1/PythonProjects/PycharmProjects/damIA/GeoDAR_ICOLD_displaced_valid/multi_pops/ICOLD2019/GeoDAR_p1975_displ.geojson'
geojson_file='/Users/ritterj1/PythonProjects/PycharmProjects/damIA/GeoDAR_ICOLD_displaced_valid/multi_pops/ICOLD2023/GeoDAR_ICOLD2023_p1975_displ.geojson'

countries_shapefile = '/Users/ritterj1/GIS_data/countries/ref-countries-2020-01m.shp/CNTR_RG_01M_2020_4326.shp/CNTR_RG_01M_2020_4326.shp'
result_gdf = get_country_iso_codes(geojson_file, countries_shapefile)
# result_gdf.loc[:,['Name of the dam','Country','ISO_CODES']] # testing


pdb.set_trace()

# Save the modified GeoDataFrame
output_file = geojson_file
result_gdf.to_file(output_file, driver='GeoJSON')