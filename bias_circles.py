import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb # pdb.set_trace()
import math

from matplotlib.cm import get_cmap
from matplotlib.collections import PathCollection
from matplotlib.legend_handler import HandlerPathCollection, HandlerLine2D

def biascircles_plot(df_bias_cntr, groupnames, grouped_by_country_counts, popgrid_shortnames, gridcolors):
    hfont = {'fontname': 'Verdana'}
    min_value = -100
    # max_value = df_bias_cntr[['gwp-un','grump-un','ghs','lscan','wpop-un']].max().max() # Venezuela has bias of 5000 % so this needs adjustment
    max_value = 1000

    scaled_angles = df_bias_cntr[['gwp-un', 'grump-un', 'ghs', 'lscan', 'wpop-un']] / max_value / 2
    scaled_angles[df_bias_cntr[['gwp-un', 'grump-un', 'ghs', 'lscan', 'wpop-un']] < 0] = -1 * df_bias_cntr[['gwp-un', 'grump-un', 'ghs', 'lscan', 'wpop-un']] / min_value / 2
    scaled_angles[df_bias_cntr[['gwp-un', 'grump-un', 'ghs', 'lscan', 'wpop-un']] > max_value] = 0.5  # Make outliers fit into plot range
    print(scaled_angles)

    r = 1.5  # outer radius of the chart
    r_inner = 0.4  # inner radius of the chart
    width = (r - r_inner) / len(popgrid_shortnames)  # width of the rings

    # Loop over countries
    for ii in range(len(df_bias_cntr)):

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.axis("equal")

        # White background circle
        circleshape0 = plt.Circle((0, 0), r, fill=True, color='white', linewidth=2, clip_on=False)
        ax.add_artist(circleshape0)

        for i, popgrid in enumerate(popgrid_shortnames):
            print(groupnames[ii], popgrid, scaled_angles.iat[ii, i])
            radius = r - i * width
            barval = scaled_angles.iat[ii, i]
            # barval = -0.05  # for legend printing only!!!
            if not math.isnan(barval):
                if barval >= 0:
                    counterclockw = False
                else:
                    counterclockw = True
                    barval = barval * (-1)
                ax.pie([barval, 1 - barval], radius=radius, startangle=90, counterclock=counterclockw, colors=[gridcolors[popgrid], 'white'],
                       wedgeprops={'width': width, 'edgecolor': 'white'})  # ,labels=[f'{cathegories[i]} – {percent[i]}%'], labeldistance=None,
            else:
                ax.text(0, radius - width / 2, 'x', ha='center', va='center', fontweight='bold', fontsize=27, **hfont)  # , fontweight='bold'

        # plt.axvline(0, color='grey', linewidth=2)#, zorder=4
        plt.plot([0, 0], [r_inner, r], color='grey', linewidth=2, clip_on=False)
        plt.plot([0, 0], [-r_inner, -r - width / 2], color='grey', linewidth=2, clip_on=False)
        circleshape1 = plt.Circle((0, 0), r_inner, fill=False, color='grey', linewidth=2, clip_on=False)
        ax.add_artist(circleshape1)
        circleshape2 = plt.Circle((0, 0), r, fill=False, color='grey', linewidth=2, clip_on=False)
        ax.add_artist(circleshape2)

        ax.text(0, 0, groupnames[ii], horizontalalignment='center', verticalalignment='center', fontweight='bold', fontsize=27, **hfont)
        ax.text(0, -r_inner / 2, str(grouped_by_country_counts[groupnames[ii]]), horizontalalignment='center', verticalalignment='center', fontweight='bold', fontsize=23, **hfont)

        plt.tight_layout()
        # plt.show()

        plt.savefig(groupnames[ii] + ".pdf", format="pdf", bbox_inches="tight", pad_inches=0.1)
        plt.savefig(groupnames[ii] + ".png", format="png", bbox_inches="tight", transparent=True, pad_inches=0.1)
