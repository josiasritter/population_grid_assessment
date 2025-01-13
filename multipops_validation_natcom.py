import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb # pdb.set_trace()
import math

from matplotlib.cm import get_cmap
from matplotlib import colors
from matplotlib.lines import Line2D

path = 'multipops_ICOLD2023_natcomm_rev.geojson'

area_adj = 'uniform' # 'none', 'uniform', 'individual'
resettle_redfact = 1    # [0.5, 1] reduces reported resettlement values
density_filt = 1500  # people/km2 to remove reservoirs with higher population densities based on ICOLD resettlement and bias-adjusted polygon area

n_grids = 5     # 5 or 8
logax = 'on'    # Switch for plotting in logarithmic scale
china = 'in'   # out or in. Switch to remove Chinese dams from analysis

gdf_in = gpd.read_file(path)
#pdb.set_trace()
if n_grids == 5:
    popgrid_names = ['WorldPop', 'LandScan', 'GHS-POP', 'GRUMP', 'GWP']
    #popgrid_names = ['LandScan', 'WorldPop 1 km', 'WorldPop 100 m', 'GHS-POP', 'GRUMP']
    popgrid_names.reverse()
    popgrid_shortnames = ['wpop-un', 'lscan', 'ghs', 'grump-un', 'gwp-un']
    #popgrid_shortnames = ['lscan', 'wpop-un1k','wpop-un', 'ghs', 'grump-un']
    popgrid_shortnames.reverse()
elif n_grids == 8:
    popgrid_names = ['GHS-POP','LandScan','WorldPop (no adjustment)','WorldPop UN-adjusted','GRUMP (no adjustment)','GRUMP UN-adjusted','GWP (no adjustment)','GWP UN-adjusted']
    popgrid_shortnames = ['ghs','lscan','wpop','wpop-un','grump','grump-un','gwp','gwp-un']

# Exclude dams with mentions on earlier year of construction (e.g. heightening etc.)
ext_dams = ['ROSEIRES', 'FOMO LOWER', 'FOMO UPPER', 'CENTIANHE', 'HUALONG', 'BUKOWKA', 'BURRINJUCK','KINCHANT','HINZE']
gdf_in = gdf_in.loc[~(gdf_in['Name of the dam'].isin(ext_dams))]

#pdb.set_trace()

# Fix area unit errors in ICOLD database (some entries use km2 but most use 10^3 km2). Divide by 1000 where ratio between values is larger than factor 100
#gdf_in = gdf_in.loc[(gdf_in['Area of Reservoir'].notnull())]
#gdf_in['Area of Reservoir'] = np.where(((gdf_in['Area of Reservoir'] / gdf_in['plg_a_km2']) > 1e4), gdf_in['Area of Reservoir'] / 1e6, gdf_in['Area of Reservoir'])   # 10^6 km2 to km2
#gdf_in['Area of Reservoir'] = np.where(((gdf_in['Area of Reservoir'] / gdf_in['plg_a_km2']) > 1e2), gdf_in['Area of Reservoir'] / 1e3, gdf_in['Area of Reservoir'])     # 10^3 km2 to km2

# Filter out dams where polygon area and ICOLD area value do not correspond
#gdf_in = gdf_in.loc[(((gdf_in['plg_a_km2']) / gdf_in['Area of Reservoir']) < 4/3) & ((gdf_in['plg_a_km2'] / (gdf_in['Area of Reservoir'])) > 2/3)]

# Exclude China to test their dominance on results due to high numbers of dams
if china == 'out':
    gdf_in = gdf_in.loc[(gdf_in['Country'] != 'China')]
    #gdf_in = gdf_in.loc[(gdf_in['Country'] != 'Brazil')]
    #dams_cou_counts = gdf_in.groupby(['Country'])['Country'].count()
    #dams_year_counts = gdf_in2.groupby(['refyear'])['refyear'].count()

#pdb.set_trace()

## Check population densities in selected reservoirs and filter out too high density reservoirs
#gdf_in['popdens_ICOLD0'] = gdf_in['Resettlement'] / gdf_in['Area of Reservoir']
#gdf_in['popdens_ICOLD0'].plot(kind='hist', bins=np.arange(0, 4600, 300))
gdf_in['popdens_ICOLD'] = gdf_in['Resettlement'] / (gdf_in['plg_a_km2']/(1-0.188))  # population density in the polygons based on ICOLD resettlement data and polygon area size accounting for polygon area bias of -18.8%
#gdf_in['popdens_ICOLD'].plot(kind='hist', bins=np.arange(0, 5000, 500))


#gdf_in = gdf_in.loc[gdf_in.popdens_ICOLD < density_filt]
gdf_in = gdf_in.loc[gdf_in.popdens_ICOLD < density_filt]
total_popdens = gdf_in['Resettlement'].sum() / (gdf_in['plg_a_km2'].sum() / (1-0.188))  #
#gdf_in['popdens_ICOLD2'].plot(kind='hist', bins=np.arange(0, 1800, 300))

#gdf_highpopdens = gdf_in.loc[gdf_in.popdens_ICOLD > 1000]   # This brings out 17 Chinese reservoirs with density > 1000 people/km2, of which 14 completed in GHS time, 3 in Grump/GHS time. Removing them would not be a big harm.
#gdf_in['popdens_ICOLD'].plot(kind='hist', bins=np.arange(0, 1100, 100))
#gdf_highpopdens = gdf_in.loc[gdf_in.popdens_ICOLD > 500]   # This brings out 53 Chinese reservoirs with density > 500 people/km2, most of them in GHS time, 1 in 2000, 1 in 2010. Removing them would not be a big harm.
#gdf_in['popdens_ICOLD'].plot(kind='hist', bins=np.arange(0, 550, 50))
#gdf_highpopdens = gdf_in.loc[gdf_in.popdens_ICOLD > 300]   # This brings out 76 Chinese and 1 Pakistani reservoir with density > 300 people/km2, most of them in GHS time, 2 in 2000, 1 in 2010. Removing them a challenge for error scores.
#gdf_in['popdens_ICOLD'].plot(kind='hist', bins=np.arange(0, 310, 10))
#gdf_highpopdens.groupby(['Country'])['Country'].count()

#pdb.set_trace()

## Bar plot of numbers of dams per map
from damnumbers import damnumbers_barplot
damnumbers_barplot(gdf_in, china)

## Area unit correction test plot
from area_scatter import area_scatterplot
#area_scatterplot(gdf_in, logax)

## Apply area adjustment factor to population estimates in order to mitigate systematic area underrepresentation in GeoDAR polygons

if area_adj == 'none':
    gdf_in['area_adj_fct'] = [1]*len(gdf_in)        # default (no adjustment)
elif area_adj == 'uniform':
    adj_fact = 0.812 / 1   # Inversed! Uniform regression-based adjustment factor (analysis of all GeoDAR dams with outlier filter 5 found a mean bias rate of -18.8%)
    gdf_in['area_adj_fct'] = 0.812 / 1  # Inversed! Uniform regression-based adjustment factor (analysis of all GeoDAR dams with outlier filter 5 found a mean bias rate of -18.8%)
elif area_adj == 'individual':
    #gdf_in['area_adj_fct'] = gdf_in['Area of Reservoir'] / gdf_in['plg_a_km2']  # Individual adjustment factor for each dam according to ratio of ICOLD reported area and polygon geometry area
    gdf_in['area_adj_fct'] = gdf_in['plg_a_km2'] / gdf_in['Area of Reservoir']   # inversed! # Individual adjustment factor for each dam according to ratio of ICOLD reported area and polygon geometry area

gdf_in['area_adj'] = gdf_in['plg_a_km2'] / 0.812

# Plot distribution of area_adj_fct
#from adj_factor_barplot import adj_factor_barplot
#adj_factor_barplot(gdf_in)

# Testing
#gdf_adj_08_12 = gdf_in.loc[(gdf_in['area_adj_fct'] > 0.8) & (gdf_in['area_adj_fct'] < 1.2)]     # This gives 125 dams in 18 countries on all continents
#gdf_adj_07_13 = gdf_in.loc[(gdf_in['area_adj_fct'] > 0.7) & (gdf_in['area_adj_fct'] < 1.3)]     # This gives 164 dams in 22 countries on all continents, might be worth doing!
#gdf_adj_07_15 = gdf_in.loc[(gdf_in['area_adj_fct'] > 0.75) & (gdf_in['area_adj_fct'] < 1.5)]     # This gives 186 dams in 27 countries on all continents. With still 28 dams in 2000-2010, compared to initially 35. This looks good!!!
#gdf_adj_lt07 = gdf_in.loc[(gdf_in['area_adj_fct'] < 0.75)]
#gdf_adj_gt15 = gdf_in.loc[(gdf_in['area_adj_fct'] > 1.5)]

# Check spatial distribution
#dams_cou_counts = gdf_s10_tight.groupby(['Country'])['Country'].count()
#dams_con_counts = gdf_s10_tight.groupby(['Continent'])['Continent'].count()

# # Testing of influence of chosen projection on polygon area calculation (all global equal-area projections)
# #gdf_adj_fct_lt1_mw = gdf_adj_fct_lt1.to_crs('ESRI:54009').copy()
# gdf_adj_fct_lt1_wcea = gdf_adj_fct_lt1.to_crs('ESRI:54034').copy() # This projection was used by GeoDAR area calculation
# #gdf_adj_fct_lt1_scea = gdf_adj_fct_lt1.to_crs('ESRI:53034').copy()
# #gdf_adj_fct_lt1_lambert = gdf_adj_fct_lt1.to_crs('IGNF:ETRS89LAEA').copy()
# # Conclusions: 1. Projection affects only by about 1-2 % (as long as any global equal-area projection is chosen). 2. GeoDAR can actually overestimate polygons (e.g. Porce II). -> Let's allow for adjustment factors < 1

#pdb.set_trace()

for popgrid in popgrid_shortnames:
    gdf_in[popgrid] = gdf_in[popgrid] / gdf_in['area_adj_fct']          # Apply bias adjustment to the dams


## Reduction factor for ICOLD Resettlement values to account for resettlement outside reservoir area
gdf_in.Resettlement = gdf_in.Resettlement.apply(lambda x: x * resettle_redfact)

#pdb.set_trace()

## Write output polygon file with results
#gdf_in.to_file('results.geojson')

### World map showing centroids of reservoirs in colour of reference year

# Read countries shapefile
world = gpd.read_file('CNTR_RG_10M_2020_4326.geojson')
world = world.loc[~(world['ISO3_CODE'].isin(['ATA']))]  # Remove Antarctica


# Colour-blind friendly scale
cmap = get_cmap('YlGn')
yearcolors = {1975: cmap(1 / 8), 1980: cmap(1.8 / 8), 1985: cmap(2.5 / 8), 1990: cmap(3.2 / 8), 1995: cmap(4 / 8), 2000: cmap(5 / 8), 2005: cmap(6.5 / 8), 2010: cmap(7.9 / 8)}
# Used before editorial requests for colour blindness
# cmap = get_cmap('gist_rainbow')
# yearcolors = {1975: cmap(0 / 8), 1980: cmap(1.1 / 8), 1985: cmap(1.5 / 8), 1990: cmap(2.1 / 8), 1995: cmap(4.1 / 8), 2000: cmap(5 / 8), 2005: cmap(6 / 8), 2010: cmap(7.5 / 8)}

from damlocations import damlocations_plot
#from damlocations_v2 import damlocations_plot
damlocations_plot(gdf_in, world, yearcolors)
pdb.set_trace()

## Prepare plot
if logax == 'on': # Increase zero values to 1, so they can be in the logarithmic plot
    gdf_in.Resettlement = gdf_in.Resettlement.apply(lambda x: 1 if x < 1 else x)
    for popgrid in popgrid_shortnames:
        gdf_in[popgrid] = gdf_in[popgrid].apply(lambda x: 1 if x < 1 else x)

if n_grids == 5:
    max_val = max(gdf_in.Resettlement.max(), gdf_in.ghs.max(), gdf_in['wpop-un'].max(), gdf_in['grump-un'].max(), gdf_in['gwp-un'].max(), gdf_in.lscan.max())
elif n_grids == 8:
    max_val = max(gdf_in.Resettlement.max(), gdf_in.ghs.max(), gdf_in.lscan.max(), gdf_in.wpop.max(), gdf_in.grump.max(), gdf_in.gwp.max(), gdf_in['wpop-un'].max(), gdf_in['grump-un'].max(), gdf_in['gwp-un'].max())
#max_val = 100000

## Performances of year 2000 against each other

fontsz = 10
fontsz2 = 11
textsp = .77

#cmap = get_cmap('Set2')
#gridcolors = {popgrid_shortnames[0]:cmap(2/8), popgrid_shortnames[1]:cmap(3/8), popgrid_shortnames[2]:cmap(4/8), popgrid_shortnames[3]:cmap(5/8), popgrid_shortnames[4]:cmap(6/8)}
cmap = get_cmap('Set1')
gridcolors = {popgrid_shortnames[0]:cmap(1/8), popgrid_shortnames[1]:cmap(7/8), popgrid_shortnames[2]:cmap(3/8), popgrid_shortnames[3]:cmap(0/8), popgrid_shortnames[4]:cmap(4/8)} #reversed
gridsyms = {popgrid_shortnames[0]:'^', popgrid_shortnames[1]:'v', popgrid_shortnames[2]:'o', popgrid_shortnames[3]:'*', popgrid_shortnames[4]:'d'}

# Plot: Performance for year 2000
from year2000_scatter import year2000_scatterplot
year2000_scatterplot(gdf_in, popgrid_names, popgrid_shortnames, gridcolors, gridsyms, logax, fontsz2)


## Plot: Performance by map year
from years_scatter import years_scatterplot
years_scatterplot(gdf_in, n_grids, popgrid_names, popgrid_shortnames, yearcolors, logax, max_val, fontsz, fontsz2)

### Plot: Performance by World Bank income class
from income_scatter import income_scatterplot
income_scatterplot(gdf_in, popgrid_names, popgrid_shortnames, logax, max_val, fontsz, fontsz2)

## Line plots of errors for each data source (and mean), seperately over map year and income class (could be additionally for different population density ranges?)

from linemetrics import linemetrics
linemetrics(gdf_in, popgrid_names, popgrid_shortnames, gridcolors, gridsyms, fontsz2)
#pdb.set_trace()

## Statistical analysis of the properties of rural areas
#pdb.set_trace()
gdf_alt = gdf_in.loc[(gdf_in['Altitude'].notnull())]
altitude = pd.to_numeric(gdf_alt['Altitude'].astype(str).str.replace(',', '.'), errors='coerce') # Change commas in original data to points
variables = [gdf_in.plg_a_km2, gdf_in.Resettlement, gdf_in.popdens_ICOLD, altitude]

def data_statistics(variables):
    # Create an empty DataFrame to store the results
    statistics_df = pd.DataFrame(columns=['Mean', 'Median', 'Std', 'Min', 'Max'])

    # Loop through each Series in the list
    for series in variables:
        # Get the name of the Series for labeling
        series_name = series.name

        # Calculate the statistical values
        mean = series.mean()
        median = series.median()
        std_dev = series.std()
        min_val = series.min()
        max_val = series.max()

        # Add the results to the DataFrame
        statistics_df.loc[series_name] = [mean, median, std_dev, min_val, max_val]

    return statistics_df

statistics_df = data_statistics(variables)
statistics_df = statistics_df.round(1)
new_row_names = ['Area [km2]', 'Population numbers [people]', 'Population density [people / km2]', 'Altitude [mamsl]']
statistics_df.index = new_row_names

statistics_df.to_csv('stat_properties.csv', header=True)#, index=False)#, mode='a')

pdb.set_trace()

## Histogram of polygon area sizes

#bins = [1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000] # Predefined bins (area ranges) for the histogram
bins = [1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000, 5000] # Predefined bins (area ranges) for the histogram
bins = [1, 5, 10, 25, 50, 100, 250, 500, 1000, 5000] # Predefined bins (area ranges) for the histogram
#import mapclassify
#bins_jenks = mapclassify.NaturalBreaks(gdf_in['plg_a_km2'], k=10)

# Call the function to plot the histogram
from polsize_barplot import polsize_barplot
#polsize_barplot(gdf_in, popgrid_shortnames, 'plg_a_km2', bins, xlabel='Surface area [km²]', ylabel='Number of rural areas evaluated')
polsize_barplot(gdf_in, popgrid_shortnames, 'area_adj', bins, xlabel='Surface area [km²]', ylabel='Number of rural areas evaluated')
pdb.set_trace()

### Circular bar charts for country-specific bias scores

# Over country
from linemetrics import error_dfs
grouped_by_country = gdf_in.groupby('ISO_CODES')
df_mape_cntr, df_smape_cntr, df_r2_cntr, df_corr_cntr, df_bias_cntr = error_dfs(grouped_by_country, popgrid_shortnames, nan_switch='off')
grouped_by_country_counts = grouped_by_country['ISO_CODES'].count()
groupnames = [name for name,unused_df in grouped_by_country]

def custom_agg(series):
    counts = series.value_counts().sort_index(ascending=False)
    result = [f"{year} ({count})" for year, count in counts.items()]
    return ', '.join(result)

# Group by "ISO_CODES" and aggregate "refyear" values with sorting
refyearstring = gdf_in.groupby('ISO_CODES')['refyear'].apply(custom_agg)
df_bias_cntr['mean'] = df_bias_cntr[popgrid_shortnames].mean(axis=1).round(decimals=1) # Calculate mean bias for each country
df_biasperc_country = df_bias_cntr.join(refyearstring)
df_biasperc_country = df_biasperc_country.reset_index()
df_biasperc_country.rename(columns = {'index':'ISO3','gwp-un':popgrid_names[0],'grump-un':popgrid_names[1],'ghs':popgrid_names[2],'lscan':popgrid_names[3],'wpop-un':popgrid_names[4],'mean':'Mean','refyear':'Reference years (number of areas evaluated)'}, inplace = True)
df_biasperc_country.to_csv('biasperc_country.csv', header=True, index=False)#, mode='a')

# PLOT
from bias_circles import biascircles_plot
biascircles_plot(df_bias_cntr, groupnames, grouped_by_country_counts, popgrid_shortnames, gridcolors)

### World Map with mean bias percentage in each country

# Merge included countries with bias dataframe
incl_countries = world[world['ISO3_CODE'].isin(groupnames)]
world_bias = incl_countries.merge(df_biasperc_country, left_on='ISO3_CODE', right_on='ISO3')    # Merge country shapes with bias scores

# Create map
from meanbias_map import meanbias_map_plot
meanbias_map_plot(world, world_bias)

pdb.set_trace()
