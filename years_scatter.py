import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb # pdb.set_trace()
import math

from matplotlib.cm import get_cmap
from scores import scores


def years_scatterplot(gdf_in, n_grids, popgrid_names, popgrid_shortnames, yearcolors, logax, max_val, fontsz, fontsz2):

    if n_grids == 5:
        rows, cols = 2, 3
        fig, ax = plt.subplots(rows, cols, figsize=(15, 10), sharey='row')  # , sharex='col')
    elif n_grids == 8:
        rows, cols = 4, 2
        fig, ax = plt.subplots(rows, cols, figsize=(11, 20), sharey='row')  # , sharex='col')

    # Loop over popgrid data sources
    for i, popgrid in enumerate(popgrid_shortnames):

        if n_grids == 5:
            if i < 2:
                ax = plt.subplot(rows, cols, i + 1)
            else:
                ax = plt.subplot(rows, cols, i + 2)  # skip last subplot in first row
        elif n_grids == 8:
            ax = plt.subplot(rows, cols, i + 1)

        # ax.set_facecolor("lightgrey")
        if logax == 'on':
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(0, 1.2 * max_val)
            ax.set_ylim(0, 1.2 * max_val)
        else:
            ax.set_xlim(0, 1.05 * max_val)
            ax.set_ylim(0, 1.05 * max_val)

        gdf_filt = gdf_in.loc[(gdf_in[popgrid_shortnames[i]].notnull())]

        # Performance scores
        # gdf_scores = gdf_filt.loc[(gdf_filt.Resettlement > 10) & (gdf_filt[popgrid_shortnames[i]] > 10)]
        MAPE, SMAPE, MAE, RMSE, R2, bias_tot, bias_mean, REE_list = scores(gdf_filt.Resettlement.values, gdf_filt[popgrid_shortnames[i]].values)

        # gdf_filt.loc[:, ['Name of the dam', 'Resettlement', 'wpop', 'grump', 'gwp', 'lscan']] # testing
        gdf_sort = gdf_filt.sort_values(by=['Resettlement'], ascending=True)
        refyears = list(np.sort(gdf_filt.refyear.unique()))

        grouped = gdf_filt.groupby('refyear')
        for key, group in grouped:
            group.plot(ax=ax, kind='scatter', x='Resettlement', s=45, y=popgrid_shortnames[i], label=key, color=yearcolors[key], zorder=3)
        # if logax == 'off':
        # ax.plot(gdf_sort.Resettlement, m * gdf_sort.Resettlement + b, color='black', linewidth=1, linestyle='dashed', alpha=0.5)    # Regression line
        ax.grid(zorder=0, which='both', color='lightgrey', alpha=0.3)  # axis='y'
        ax.plot([0, 1.2 * max_val], [0, 1.2 * max_val], color='grey', linewidth=1, linestyle='dashed', alpha=0.8, zorder=0)  # 1-1 line
        ax.text(.04 * max_val, 1, '$N$ = ' + str(len(gdf_filt)) + '\n$Bias$ = ' + str(bias_tot) + '%' + '\n$sMAPE$ = ' + str(SMAPE), ha='left', va='bottom', color='black', fontsize=fontsz)  # , transform=ax.transAxes)
        # ax.text(.04 * max_val, 1, '$N$ = '+str(len(gdf_filt)) + '\n$R^2$ = '+str(R2) + '\n$MAPE$ = '+str(MAPE) + '\nMean bias = '+str(bias_mean) + '\nTotal bias = '+str(bias_tot), ha='left', va='bottom', color='black', fontsize=fontsz)#, transform=ax.transAxes)
        ax.set_xlabel("Reported population $P_{reported}$", labelpad=1, fontsize=fontsz2)  # , fontsize=fontsz
        if n_grids == 5 and i in [0, 2]:
            ax.set_ylabel("Predicted population $P_{predicted}$", labelpad=1, fontsize=fontsz2)  # , fontsize=fontsz
        elif n_grids == 8 and i in [0, 2, 4, 6]:
            ax.set_ylabel("Predicted population $P_{predicted}$", labelpad=1, fontsize=fontsz2)  # , fontsize=fontsz
        ax.set_title(popgrid_names[i], fontweight='demibold')

        ax.legend(title="Map year")
        # lgnd0 = ax.legend(loc='upper left')
        # lgnd0.legendHandles[0]._sizes = [40]
        # lgnd0.legendHandles[1]._sizes = [40]
        # lgnd0.legendHandles[2]._sizes = [40]

    plt.tight_layout()
    # plt.savefig("allyears_uniformadj.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    # pdb.set_trace()
