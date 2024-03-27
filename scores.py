import geopandas as gpd
import pandas as pd
import numpy as np
import pdb # pdb.set_trace()

# Performance scores
from sklearn import metrics    # (https://scikit-learn.org/stable/modules/model_evaluation.html)
def scores(rep, pred):
    MAPE = np.around(metrics.mean_absolute_percentage_error(rep[rep != 0], pred[rep != 0]), decimals=2)
    SMAPE = np.around(np.mean(np.abs(pred - rep) / ((np.abs(pred) + np.abs(rep)))), decimals=2) # simple/unscaled SMAPE (range 0-1)
    #SMAPE = np.around(np.mean(np.abs(pred - rep) / ((np.abs(pred) + np.abs(rep)) / 2)), decimals=2) # scaled SMAPE (range 0-2)
    bias_mean = np.around(np.mean(pred[rep != 0] / rep[rep != 0]), decimals=2)
    MAE = np.around(metrics.mean_absolute_error(rep, pred), decimals=2)
    RMSE = np.around(metrics.mean_squared_error(rep, pred, squared=False), decimals=2)
    R2 = np.around(metrics.r2_score(rep, pred), decimals=2)
    #bias_tot = np.around(np.sum(pred)/np.sum(rep), decimals=2)
    bias_tot = np.around(100*((np.sum(pred) - np.sum(rep)) / np.sum(rep)), decimals=1)
    REE_list = (pred - rep) / rep * 100  # [%]
    return MAPE, SMAPE, MAE, RMSE, R2, bias_tot, bias_mean, REE_list