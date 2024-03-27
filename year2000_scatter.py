import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb # pdb.set_trace()
import math

from matplotlib.cm import get_cmap

def year2000_scatterplot(gdf_in, popgrid_names, popgrid_shortnames, gridcolors, gridsyms, logax, fontsz2):
    fig, ax = plt.subplots(1, 1, figsize=(13, 5))

    gdf_filt = gdf_in.loc[(gdf_in.refyear == 2000)]
    # MAPE, SMAPE, MAE, RMSE, R2, bias_tot, bias_mean, REE_list = scores(gdf_filt.Resettlement.values, gdf_filt[popgrid_shortnames[i]].values) # Would stats make sense here??? Some yes, e.g. total bias, mean bias,
    # grouped = gdf_filt.groupby('wb_inclvl2000')
    max_val2000 = max(gdf_filt.Resettlement.max(), gdf_filt.ghs.max(), gdf_filt.wpop.max(), gdf_filt.grump.max(), gdf_filt.gwp.max(), gdf_filt.lscan.max())

    if logax == 'on':
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0, 1.2 * max_val2000)
        ax.set_ylim(0, 1.2 * max_val2000)

    for i, popgrid in enumerate(popgrid_shortnames):
        gdf_filt.plot(ax=ax, kind='scatter', x='Resettlement', s=45, y=popgrid, label=popgrid_names[i], color=gridcolors[popgrid], marker=gridsyms[popgrid], zorder=3)

    ax.grid(zorder=0, which='both', color='lightgrey', alpha=0.3)  # axis='y'
    ax.plot([0, 1.2 * max_val2000], [0, 1.2 * max_val2000], color='grey', linewidth=1, linestyle='dashed', alpha=0.8, zorder=0)  # 1-1 line
    # ax.text(textsp * max_val2000, 0.01, '$N$ = '+str(len(rep_area)) + '\n$MAPE$ = '+str(MAPE_area) + '\n$R^2$ = '+str(R2_area) + '\n$Mean$ $bias$ = '+str(bias_mean_area) + '\n$Total$ $bias$ = '+str(bias_tot_area), ha='left', va='bottom', color='black', fontsize=fontsz)#, transform=ax.transAxes)
    # ax.text(.11 * max_val2000, 1, '$N$ = ' + str(len(gdf_filt)) + '\n$R^2$ = ' + str(R2) + '\nBias = ' + str(bias_tot), ha='left', va='bottom', color='black', fontsize=fontsz)  # , transform=ax.transAxes)
    # ax.set_xlabel("Reported resettled people $R_{reported}$", labelpad=1, fontsize=fontsz2)  # , fontsize=fontsz
    # ax.set_ylabel("Predicted resettled people $R_{predicted}$", labelpad=1, fontsize=fontsz2)  # , fontsize=fontsz
    ax.set_xlabel("Reported population $P_{reported}$", labelpad=1, fontsize=fontsz2)  # , fontsize=fontsz
    ax.set_ylabel("Predicted population $P_{predicted}$", labelpad=1, fontsize=fontsz2)  # , fontsize=fontsz
    ax.set_title('Year 2000', fontweight='demibold')

    ax.legend(title="Population map", loc='upper left')

    plt.tight_layout()
    plt.show()

    # pdb.set_trace()