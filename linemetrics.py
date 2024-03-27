import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb # pdb.set_trace()
import math

from matplotlib.cm import get_cmap
from scores import scores

def error_dfs(grouped, popgrid_shortnames, nan_switch='on'):  # generates error summary dataframes
    # Create empty dataframes to be filled with error metrics
    groupnames = [name for name, unused_df in grouped]
    df_mape = pd.DataFrame(np.nan, index=groupnames, columns=popgrid_shortnames)
    df_smape = df_mape.copy()
    df_bias = df_mape.copy()
    df_r2 = df_mape.copy()
    df_corr = df_mape.copy()

    # Calculate scores and fill into dataframes
    for key, group in grouped:
        for i, popgrid in enumerate(popgrid_shortnames):
            if nan_switch == 'on':
                if not group[popgrid].isnull().any().any():
                    MAPE, SMAPE, MAE, RMSE, R2, bias_tot, bias_mean, REE_list = scores(group.Resettlement.values, group[popgrid].values)
                    df_mape.at[key, popgrid] = MAPE  # Fill value into dataframe
                    df_smape.at[key, popgrid] = SMAPE
                    df_r2.at[key, popgrid] = R2
                    df_bias.at[key, popgrid] = bias_tot
                    df_corr.at[key, popgrid] = group['Resettlement'].corr(group[popgrid], method='spearman')
            if nan_switch == 'off':  # Works only for bias_percentage
                group_notnan = group[~group[popgrid].isna()]
                bias_tot = np.around(100 * ((group_notnan[popgrid].sum() - group_notnan['Resettlement'].sum()) / group_notnan['Resettlement'].sum()), decimals=1)
                df_bias.at[key, popgrid] = bias_tot

    if nan_switch == 'on':
        # Add mean
        df_mape['mean'] = df_mape.mean(axis=1, skipna=True)  # Based on MAPE, GHS looks much better than it is, simply because negative errors are limited by 100%, but positive errors can go to inf.
        df_smape['mean'] = df_smape.mean(axis=1, skipna=True)
        df_r2['mean'] = df_r2.mean(axis=1, skipna=True)
        df_bias['mean'] = df_bias.mean(axis=1, skipna=True)
        df_corr['mean'] = df_corr.mean(axis=1, skipna=True)

    return df_mape, df_smape, df_r2, df_corr, df_bias

def linemetrics(gdf_in, popgrid_names, popgrid_shortnames, gridcolors, gridsyms, fontsz2):
    # Remove values with Resettlement < 1 to avoid artifacts for SMAPE
    gdf_filt = gdf_in.loc[(gdf_in.Resettlement > 1)]

    # gdf_filt = gdf_in.copy()

    # Over map year
    grouped_by_year = gdf_filt.groupby('refyear')
    df_mape_year, df_smape_year, df_r2_year, df_corr_year, df_bias_year = error_dfs(grouped_by_year, popgrid_shortnames)
    groupnames = [name for name, unused_df in grouped_by_year]
    # for df_it in [df_mape_year, df_r2_year, df_bias_year]:      # Remove mean where only one data source
    #    df_it.at[groupnames[0:2], 'mean'] = np.nan
    #    df_mape_year.at[groupnames[0:2], 'mean'] = NaN

    # Over country income class
    gdf_filt = gdf_filt.loc[gdf_filt.refyear >= 2000]  # remove entries before year 2000 to avoid dominance of GHS-POP and GRUMP
    # grouped_by_income = gdf_filt.groupby('wb_inclvl2000')
    grouped_by_income = gdf_filt.groupby('wb_inclvl_refyear')
    df_mape_income, df_smape_income, df_r2_income, df_corr_income, df_bias_income = error_dfs(grouped_by_income, popgrid_shortnames)

    # Line plots
    # rows, cols = 2, 2
    # fig, ax = plt.subplots(rows, cols, figsize=(11, 6), sharey='row', sharex='col')
    rows, cols = 2, 1
    fig, ax = plt.subplots(rows, cols, figsize=(6, 6), sharey='row', sharex='col')
    linecols = [*gridcolors.values()] + ['black']
    linestyles = [i + '-' for i in [*gridsyms.values()]] + ['s--']

    # Loop over figure panels
    # for i, error_df in enumerate([df_bias_year, df_bias_income, df_smape_year, df_smape_income]): # 4 panels
    for i, error_df in enumerate([df_bias_year, df_smape_year]):  # only year
        # for i, error_df in enumerate([df_bias_income, df_smape_income]):    # only income level
        ax = plt.subplot(rows, cols, i + 1)
        error_df.plot(ax=ax, color=linecols, style=linestyles)
        ax.grid(axis='y', zorder=0, color='lightgrey', alpha=0.3)
        if i == 0:
            ax.set_ylabel("Bias percentage", labelpad=1, fontsize=fontsz2)
            handles, labels = plt.gca().get_legend_handles_labels()
            ax.legend(handles, popgrid_names + ['Mean of all'], loc='upper left')  # , title="Data source" # 4 panels or only year
            # ax.legend([handles[0]] + handles[2:], [popgrid_names[0]] + popgrid_names[2:] + ['Mean of all'], loc='lower left') #, loc='upper right', bbox_to_anchor=(0,-88), title="Data source"    # only income level
        # if i == 2: # 4 panels
        if i == 1:  # only year or only income level
            ax.set_xlabel("Map year", labelpad=2, fontsize=fontsz2)  # 4 panels or only year
            # ax.set_xlabel("Country income level", labelpad=2, fontsize=fontsz2) # only income level
            ax.set_ylabel("$sMAPE$", labelpad=1, fontsize=fontsz2)
        if i == 3:
            ax.set_xlabel("Country income level", labelpad=2, fontsize=fontsz2)
        if i != 0:
            ax.legend('', frameon=False)

    plt.tight_layout()
    plt.show()
