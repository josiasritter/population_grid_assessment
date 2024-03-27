import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb # pdb.set_trace()
import math

from matplotlib.cm import get_cmap

def damnumbers_barplot(gdf_in, china):
    if china == 'in':
        width = 0.5
        ax = gdf_in.refyear.value_counts(sort=False).sort_index().plot(kind='bar', color='dimgrey', width=width, alpha=1, label='All countries')  # All dams
        ax.bar_label(ax.containers[0], fontsize=8, fontweight='bold', color='dimgrey', alpha=1)

        gdf_china = gdf_in.loc[(gdf_in['Country'] == 'China')]
        ax2 = gdf_china.refyear.value_counts(sort=False).sort_index().plot(kind='bar', color='silver', width=0.7 * width, alpha=1, label='China')  # China
        # ax2.bar_label(ax2.containers[0], fontsize=8, fontweight='bold', color='lightgrey', alpha=1)

        plt.xlabel("Reference year of population map", fontsize=10, labelpad=5)
        # plt.ylabel("Number of reservoirs", fontsize=10, labelpad=5)
        plt.ylabel("Number of rural areas evaluated", fontsize=10, labelpad=5)

        plt.legend()

        plt.tight_layout()
        plt.show()

    # pdb.set_trace()