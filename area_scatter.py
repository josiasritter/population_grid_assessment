import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb # pdb.set_trace()
import math

from matplotlib.cm import get_cmap

def area_scatterplot(gdf_in, logax):
    max_val = 1.05 * max(max(gdf_in['plg_a_km2']), max(gdf_in['Area of Reservoir']))

    from matplotlib.cm import get_cmap
    cmap = get_cmap('gist_rainbow')
    srccolors = {'GRanD v1.3': cmap(6 / 8), 'HydroLAKES v1.0': cmap(8 / 8), 'UCLA Circa 2015': cmap(3 / 8)}

    f, ax = plt.subplots()
    # plt.scatter(valid_data, sim_data,  color='black')

    if logax == 'on':
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim(0, 1.2 * max_val)
        plt.ylim(0, 1.2 * max_val)
    else:
        plt.xlim(0, 1.05 * max_val)
        plt.ylim(0, 1.05 * max_val)

    grouped = gdf_in.groupby('plg_src')
    for key, group in grouped:
        group.plot(ax=ax, kind='scatter', x='Area of Reservoir', y='plg_a_km2', s=2, label=key, color=srccolors[key])

    # plt.plot(list(h_rep), list(WL_pred), color='crimson', linewidth=2)
    plt.plot([0, 1.2 * max_val], [0, 1.2 * max_val], color='lightgrey', linewidth=1, linestyle='dashed', alpha=0.5)
    # plt.text(.1*max_val, 1, '$N$ = '+str(len(valid_data)) + '\n$MAPE$ = '+str(MAPE)+ '\n$R^2$ = '+str(R2) + '\n$RMSE$ = '+str(RMSE) + ' $km^2$' + '\nMean bias = '+str(bias_mean) + '\nTotal bias = '+str(bias_tot), ha='left', va='bottom', color='black', fontsize=10)#, transform=ax.transAxes)
    # plt.text(.05*max_val, 1, '$N$ = '+str(len(valid_data)) + '\n$R^2$ = '+str(R2) + '\n$MAPE$ = '+str(MAPE) + '\nMean bias = '+str(bias_mean) + '\nTotal bias = '+str(bias_tot), ha='left', va='bottom', color='black', fontsize=10)#, transform=ax.transAxes)
    plt.grid(zorder=0, which='both', color='lightgrey', alpha=0.3)  # axis='y'
    plt.xlabel("Reservoir area reported by ICOLD  [$km^2$]", fontsize=10, labelpad=1)
    plt.ylabel("Reservoir area delineated from polygon [$km^2$]", fontsize=10, labelpad=1)
    plt.legend(title="Reservoir polygon source")

    plt.tight_layout()
    plt.show()

    #pdb.set_trace()
