import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb # pdb.set_trace()
import math

from matplotlib.cm import get_cmap
from matplotlib.collections import PathCollection
from matplotlib.legend_handler import HandlerPathCollection, HandlerLine2D

def damlocations_plot(gdf_in, gdf_countries, yearcolors):

    #path = '/Users/ritterj1/PythonProjects/PycharmProjects/damIA/GeoDAR_ICOLD_displaced_valid/multi_pops/ICOLD2023/multipops_ICOLD2023.geojson' # this has already filters applied: Built after 1975, Resettlement data, polygon area > 1km2, no transb.
    #gdf_in = gpd.read_file(path)

    ### World map showing centroids of reservoirs in colour of reference year

    # Compute reservoir centroids
    gdf_in_proj = gdf_in.to_crs('ESRI:54030')
    centroids = gdf_in_proj.copy()
    centroids.geometry = gdf_in_proj.centroid

    # Plot the world in light grey without shape outline
    fig, ax = plt.subplots(figsize=(13, 8))
    gdf_countries.to_crs('ESRI:54030').plot(ax=ax, color='lightgrey', edgecolor='white', linewidth=.2)

    # Overplot reservoir centroids. Markers coloured by reference year. Marker size in function of number of resettled people
    plt.rcParams['font.family'] = 'Verdana'
    centroids['colour'] = [yearcolors[i] for i in centroids['refyear']] # add color column
    centroids.Resettlement = centroids.Resettlement.apply(lambda x: 1 if x < 1 else x) # Increase Resettlement values of 0 to 1, so they can be in logarithmic scale
    centroids['size'] = [(3+0.00002*i+math.log(i)) for i in centroids['Resettlement']] # add markersize column
    #centroids['size'] = [(3+0.00001*i+math.log2(i)/math.log2(1.5))/2 for i in centroids['Resettlement']]

    #centroids.plot(ax=ax, color=centroids['colour'], markersize=centroids['size'], edgecolor="black", linewidth=0.2, legend=True)#, legend_kwdsdict={})

    groupplots = []
    grouped = centroids.groupby('refyear')
    for key, group in grouped:
        p = group.plot(ax=ax, label=key, color=yearcolors[key], marker='o', markersize=group['size'].values, edgecolor='black', linewidth=0.2) # , x='Resettlement', s=45, y=popgrid_shortnames[i], zorder=3
        groupplots.append(p)
    #pdb.set_trace()

    #lgnd = ax.legend(title="Reference year")#, markersize=(math.log(centroids['Resettlement'].max())*2)

    #lgnd.legendHandles[0]._sizes = [30]
    #lgnd.legendHandles[1]._sizes = [30]
    #lgnd.legendHandles[0]._legmarker.set_markersize(30)
    #handles, labels = ax.get_legend_handles_labels()

    def updatescatter(handle, orig):
        handle.update_from(orig)
        handle.set_sizes([30])

    def updateline(handle, orig):
        handle.update_from(orig)
        handle.set_markersize(8)

    #lgnd1 = ax.legend(title="Year", loc='upper right', handler_map={PathCollection: HandlerPathCollection(update_func=updatescatter), plt.Line2D: HandlerLine2D(update_func=updateline)}, title_fontsize=9, alignment='left', fontsize=8, handletextpad=0.1, frameon=False)
    lgnd1 = ax.legend(title="Reference year", bbox_to_anchor=(0.87, 0.8), loc='upper left', handler_map={PathCollection : HandlerPathCollection(update_func=updatescatter), plt.Line2D : HandlerLine2D(update_func = updateline)}, title_fontsize=9, alignment='left', fontsize=8, handletextpad=0.1, frameon=False)
    #for count, legend_handle in enumerate(handles):
    #    legend_handle.set(markersize = 10)#, alpha = 0.8)
    #for handle in handles:
    #    handle.set_sizes([10])
    ax.add_artist(lgnd1)

    # Add second legend related to marker size
    from matplotlib.lines import Line2D
    #labelnrs = [1,10,100,1000,10000,100000,1000000]
    #labels = [str(i) for i in labelnrs]
    #markersizes = [5+math.log(i)*2 for i in labelnrs]
    scalefac = 0.5
    legend_elements = [Line2D([0], [0], marker='o', color='black', label='100', markerfacecolor='black', markersize=math.log(10), ls = '', linewidth=0.2),
                       Line2D([0], [0], marker='o', color='black', label='10,000', markerfacecolor='black', markersize=math.log(30), ls = '', linewidth=0.2),
                       Line2D([0], [0], marker='o', color='black', label='1,000,000', markerfacecolor='black', markersize=math.log(100), ls = '', linewidth=0.2)]
    #legend_elements = []
    #for ii, label in enumerate(labels):
    #    legend_elements.append(Line2D([0], [0], marker='o', color='black', label=label, markerfacecolor='none', markersize=markersizes[ii], ls = ''))

    #lgnd2 = ax.legend(title='Reported population', handles=legend_elements, bbox_to_anchor =(0.6,0.04), loc='lower center', ncol=len(legend_elements), title_fontsize=9, fontsize=8.5, frameon=False)
    lgnd2 = ax.legend(title='Reported population', handles=legend_elements, bbox_to_anchor=(0.87, 0.97), loc='upper left', title_fontsize=9, alignment='left', fontsize=8, handletextpad=0.1, frameon=False)

    ax.set_axis_off()
    plt.show()
    #plt.savefig('/Users/ritterj1/PythonProjects/PycharmProjects/damIA/GeoDAR_ICOLD_displaced_valid/multi_pops/ICOLD2023/lowfilters/damlocations2.pdf', bbox_inches='tight', pad_inches=0.1)
