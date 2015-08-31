from parse_data import df, nbd_residues
import numpy as np
from copy import deepcopy
from matplotlib import pyplot as plt
from tbidbaxlipo.util.calculate_error_variance import calc_err_var

# There are a number of outliers in the FRET data that Justin attributes
# to debris floating in the cuvette that produces momentary spikes. Here
# we explicitly remove several of these from the 54C data, setting them to
# NaN.
df_pre = deepcopy(df)

# Outliers in 54C FRET data
outliers_54_FRET = [21, 24, 77, 103, 108, 113]
for outlier_ix in outliers_54_FRET:
    df_pre[('Bid', 'FRET', '54', 1, 'VALUE')][outlier_ix] = np.nan

# Outliers in 54C NBD data
outliers_54_NBD = [109]
for outlier_ix in outliers_54_NBD:
    df_pre[('Bid', 'NBD', '54', 1, 'VALUE')][outlier_ix] = np.nan

# NOTE: There are currently no really major outliers for 126C.

# Matrices to use in fitting
time_54 = df_pre[('Bid', 'FRET', '54', 1, 'TIME')].values
data_54 = np.zeros((1, 2, len(time_54)))
data_54[0, 0, :] = df_pre[('Bid', 'NBD', '54', 1, 'VALUE')].values
data_54[0, 1, :] = df_pre[('Bid', 'FRET', '54', 1, 'VALUE')].values

time_126 = df_pre[('Bid', 'FRET', '126', 1, 'TIME')].values
data_126 = np.zeros((1, 2, len(time_126)))
data_126[0, 0, :] = df_pre[('Bid', 'NBD', '126', 1, 'VALUE')].values
data_126[0, 1, :] = df_pre[('Bid', 'FRET', '126', 1, 'VALUE')].values

# Calculate the prior on standard error of the data by running polynomial
# fits on the final points
# 54
data_54_sigma = np.zeros((1, 2))
(residuals, _) = calc_err_var(data_54[0, 0, :], last_n_pts=80,
                              fit_type='cubic', plot=False)
data_54_sigma[0, 0] = np.std(residuals, ddof=1)

(residuals, _) = calc_err_var(data_54[0, 1, :], last_n_pts=80,
                              fit_type='cubic', plot=False)
data_54_sigma[0, 1] = np.std(residuals, ddof=1)

# 126
data_126_sigma = np.zeros((1, 2))
(residuals, _) = calc_err_var(data_126[0, 0, :], last_n_pts=80,
                              fit_type='cubic', plot=False)
data_126_sigma[0, 0] = np.std(residuals, ddof=1)

(residuals, _) = calc_err_var(data_126[0, 1, :], last_n_pts=80,
                              fit_type='cubic', plot=False)
data_126_sigma[0, 1] = np.std(residuals, ddof=1)


def plot_outliers():
    plt.ion()
    # 54C FRET
    plt.figure()
    plt.plot(df[('Bid', 'FRET', '54', 1, 'VALUE')].values, 'r')
    plt.plot(df_pre[('Bid', 'FRET', '54', 1, 'VALUE')].values, 'b')
    plt.title('54C FRET')

    # 54C NBD
    plt.figure()
    plt.plot(df[('Bid', 'NBD', '54', 1, 'VALUE')].values, 'r')
    plt.plot(df_pre[('Bid', 'NBD', '54', 1, 'VALUE')].values, 'b')
    plt.title('54C NBD')

    # 126C FRET
    plt.figure()
    plt.plot(df[('Bid', 'FRET', '126', 1, 'VALUE')].values, 'r')
    plt.plot(df_pre[('Bid', 'FRET', '126', 1, 'VALUE')].values, 'b')
    plt.title('126C FRET')

    # 126C NBD
    plt.figure()
    plt.plot(df[('Bid', 'NBD', '126', 1, 'VALUE')].values, 'r')
    plt.plot(df_pre[('Bid', 'NBD', '126', 1, 'VALUE')].values, 'b')
    plt.title('126C NBD')

if __name__ == '__main__':
    plot_outliers()
