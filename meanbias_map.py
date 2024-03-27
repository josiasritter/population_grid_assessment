import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb # pdb.set_trace()
import math

from matplotlib.cm import get_cmap
from matplotlib.collections import PathCollection
from matplotlib.legend_handler import HandlerPathCollection, HandlerLine2D

def meanbias_map_plot(world, world_bias):

    # Plot the world in light grey without shape outline
    fig, ax = plt.subplots(figsize=(13, 8))
    plt.rcParams['font.family'] = 'Verdana'
    world.to_crs('ESRI:54030').plot(ax=ax, color='lightgrey', edgecolor='none')#, edgecolor='darkgrey', linewidth=.15)  # )

    ## VERSION 1 WITH COLORBAR
    # from matplotlib.colors import LinearSegmentedColormap
    # from matplotlib.cm import ScalarMappable
    # cvals  = [-100., 0, 100.]
    # norm=plt.Normalize(min(cvals),max(cvals))
    # cmap0=get_cmap('viridis')
    # colors = [cmap0(0), cmap0(0.8), cmap0(0.99)]  # R -> G -> B
    # #n_bins = 1000  # Discretizes the interpolation into bins
    # tuples = list(zip(map(norm,cvals), colors))
    # cmap_bias = LinearSegmentedColormap.from_list("", tuples)

    # Plot the included countries coloured by mean bias
    # world_bias.to_crs('ESRI:54030').plot(ax=ax, column='Mean', cmap=cmap_bias, edgecolor='white', linewidth=.1, norm=norm)#, legend=True, font='Verdana'

    # Plot legend
    # #cbar_ax = fig.add_axes([0.85, 0.28, 0.01, 0.44])
    # cbar_ax = fig.add_axes([0.45, 0.28, 0.25, 0.015])  # horizontal alignment
    # cb = plt.colorbar(ScalarMappable(norm=norm, cmap=cmap_bias), cax=cbar_ax, orientation='horizontal')
    # cb.set_label('Mean bias percentage', fontsize=9)
    # cb.ax.tick_params(labelsize=8.5)

    ## VERSION 2 WITH DISCRETE COLORS

    # Add colour column
    #bins = [-100, -50, -25, -10, 10, 25, 50, np.inf] # Classes as in Kuffer
    #labels = ['Great underestimation\n$Bias_{mean}$ < -50%', 'Underestimation\n-50% < $Bias_{mean}$ < -25%', 'Slight underestimation\n-25% < $Bias_{mean}$ < -10%', 'Accurate\n-10% < $Bias_{mean}$ < 10%', 'Slight overestimation\n10% < $Bias_{mean}$ < 25%','Overestimation\n25% < $Bias_{mean}$ < 50%', 'Great overestimation\n50% < $Bias_{mean}$']
    # cmap = plt.get_cmap('tab20c')
    # meanbias_cols =[cmap(4/20),cmap(5/20),cmap(6/20),cmap(9/20),cmap(14/20),cmap(13/20),cmap(12/20)]
    # cmap = plt.get_cmap('tab20b')
    # meanbias_cols =[cmap(12/20),cmap(13/20),cmap(14/20),cmap(4/20),cmap(2/20),cmap(1/20),cmap(0/20)]
    cmap_r = plt.get_cmap('Reds')
    cmap_g = plt.get_cmap('Greens')
    cmap_y = plt.get_cmap('viridis')
    cmap_b = plt.get_cmap('Blues')
    #meanbias_cols = [cmap_r(.8), cmap_r(.65), cmap_r(.5), cmap_g(.5), cmap_b(.5), cmap_b(.65), cmap_b(.8)]

    bins = [-100, -50, -25, 25, 50, np.inf]     # Reduced classes as in Bai
    labels = ['Great underestimation\n$Bias_{mean}$ < -50%', 'Underestimation\n-50% < $Bias_{mean}$ < -25%', 'Accurate\n-25% < $Bias_{mean}$ < 25%',
              'Overestimation\n25% < $Bias_{mean}$ < 50%', 'Great overestimation\n50% < $Bias_{mean}$']
    #meanbias_cols = [cmap_r(.77), cmap_r(.55), cmap_g(.6), cmap_b(.55), cmap_b(.77)]
    meanbias_cols = [cmap_r(.77), cmap_r(.55), cmap_y(.99), cmap_b(.55), cmap_b(.77)]

    # Create a new column 'Color' based on the bins
    world_bias['Color'] = pd.cut(world_bias['Mean'], bins, labels=meanbias_cols, right=False)

    # Plot the included countries coloured by mean bias
    world_bias.to_crs('ESRI:54030').plot(ax=ax, color=world_bias['Color'], edgecolor='white', linewidth=.2)#, alpha=0.8)  # , legend=True, font='Verdana' # VERSION 1

    # Plot legend # VERSION 2
    import matplotlib.patches as mpatches
    #legend_elements = [mpatches.Patch(color=meanbias_cols[0], label=labels[0]),mpatches.Patch(color=meanbias_cols[1], label=labels[1]),mpatches.Patch(color=meanbias_cols[2], label=labels[2]),mpatches.Patch(color=meanbias_cols[3], label=labels[3]),mpatches.Patch(color=meanbias_cols[4], label=labels[4]),mpatches.Patch(color=meanbias_cols[5], label=labels[5]),mpatches.Patch(color=meanbias_cols[6], label=labels[6])]
    legend_elements = []
    for i, col in enumerate(meanbias_cols):
       legend_elements.append(mpatches.Patch(color=col, label=labels[i]))
    ax.legend(handles=legend_elements, bbox_to_anchor=(0.54, -0.04), loc='lower center', ncol=len(legend_elements), columnspacing=0.8, fontsize=8, handletextpad=0.5, handlelength=1, frameon=False)

    ## Finish plot

    ax.set_axis_off()
    # plt.tight_layout()
    plt.savefig('/Users/ritterj1/PythonProjects/PycharmProjects/damIA/GeoDAR_ICOLD_displaced_valid/multi_pops/ICOLD2023/lowfilters/circlebias/0_biasmap3_raw.pdf', bbox_inches='tight', pad_inches=0.1)

    #plt.show()