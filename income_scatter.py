import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb # pdb.set_trace()
import math

from matplotlib.cm import get_cmap
from scores import scores

def income_scatterplot(gdf_in, popgrid_names, popgrid_shortnames, logax, max_val, fontsz, fontsz2):
    # Rename income level abbreviations
    income_remap = {'L': 'Low', 'LM': 'Lower Middle', 'UM': 'Upper Middle', 'H': 'High'}
    for key in income_remap:
        # gdf_in['wb_inclvl2000'] = gdf_in['wb_inclvl2000'].apply(lambda x: income_remap[key] if x == key else x)
        gdf_in['wb_inclvl_refyear'] = gdf_in['wb_inclvl_refyear'].apply(lambda x: income_remap[key] if x == key else x)

    cmap = get_cmap('coolwarm')
    incomecolors = {'Low': cmap(0 / 3), 'Lower Middle': cmap(1 / 3), 'Upper Middle': cmap(2 / 3), 'High': cmap(3 / 3)}

    rows, cols = 2, 3
    fig, ax = plt.subplots(rows, cols, figsize=(15, 10), sharey='row')  # , sharex='col')

    # Loop over popgrid data sources
    for i, popgrid in enumerate(popgrid_shortnames):
        print(popgrid)
        if i < 2:
            ax = plt.subplot(rows, cols, i + 1)
        else:
            ax = plt.subplot(rows, cols, i + 2)  # skip last subplot in first row

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
        MAPE, SMAPE, MAE, RMSE, R2, bias_tot, bias_mean, REE_list = scores(gdf_filt.Resettlement.values, gdf_filt[popgrid_shortnames[i]].values)

        # grouped = gdf_filt.groupby('wb_inclvl2000')
        grouped = gdf_filt.groupby('wb_inclvl_refyear')
        for key, group in grouped:
            group.plot(ax=ax, kind='scatter', x='Resettlement', s=45, y=popgrid_shortnames[i], label=key, color=incomecolors[key], zorder=3)
        # if logax == 'off':
        #    ax.plot(gdf_sort.Resettlement, m * gdf_sort.Resettlement + b, color='black', linewidth=1, linestyle='dashed', alpha=0.5)    # Regression line
        ax.grid(zorder=0, which='both', color='lightgrey', alpha=0.3)  # axis='y'
        ax.plot([0, 1.2 * max_val], [0, 1.2 * max_val], color='grey', linewidth=1, linestyle='dashed', alpha=0.8, zorder=0)  # 1-1 line
        ax.text(.04 * max_val, 1, '$N$ = ' + str(len(gdf_filt)) + '\n$Bias$ = ' + str(bias_tot) + '%' + '\n$sMAPE$ = ' + str(SMAPE), ha='left', va='bottom', color='black', fontsize=fontsz)  # , transform=ax.transAxes)
        ax.set_xlabel("Reported population $P_{reported}$", labelpad=1, fontsize=fontsz2)  # , fontsize=fontsz
        if i in [0, 2]:
            ax.set_ylabel("Predicted population $P_{predicted}$", labelpad=1, fontsize=fontsz2)  # , fontsize=fontsz
        ax.set_title(popgrid_names[i], fontweight='demibold')

        # reordering the labels
        handles, labels = plt.gca().get_legend_handles_labels()
        if len(handles) == 4:
            order = [1, 2, 3, 0]
            ax.legend([handles[j] for j in order], [labels[k] for k in order], title="Country income level")
        else:
            ax.legend(title="Country income level")

    plt.tight_layout()
    # plt.savefig("incomeclass_uniformadj.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    # pdb.set_trace()
