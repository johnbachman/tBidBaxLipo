from tbidbaxlipo.util.calculate_error_variance import calc_err_var
import numpy as np


def mod_zscore(residuals):
    numer = 0.6745 * (residuals - np.median(residuals))
    denom = np.median(np.abs(residuals - np.median(residuals)))
    return numer / denom

def filtered(x, zscores, cutoff):
    x_copy = np.copy(x)
    for i in range(0, len(zscores)):
        if np.abs(zscores[-i]) > cutoff:
            x_copy[-i] = np.nan
    return x_copy


def find_outliers(values):
    (res, fig) = calc_err_var(values, last_n_pts=110, fit_type='cubic',
                              plot=False)
    zscores = mod_zscore(res)
    return filtered(values, zscores, 3.5)

#from tbidbaxlipo.plots.bid_bim_fret_nbd_release.preprocess_data import df
#from matplotlib import pyplot as plt
#time = df[('Bim', 'FRET', '3', 1, 'TIME')].values
#values = df[('Bim', 'FRET', '3', 1, 'VALUE')].values
#plt.ion()
#plt.figure()
#plt.plot(time, values, linestyle='', marker='.', color='r')
#plt.plot(time, find_outliers(values), linestyle='', marker='.',
#         color='b')


