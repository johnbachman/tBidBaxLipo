from parse_data import df
import numpy as np

# There are a number of outliers in the FRET data that Justin attributes
# to debris floating in the cuvette that produces momentary spikes. Here
# we explicitly remove several of these from the 54C data, setting them to
# NaN.
outliers = [21, 24, 77, 103, 108, 113]
for outlier_ix in outliers:
    df[('54', 'FRET', 'VALUE')][outlier_ix] = np.nan

# Matrices to use in fitting

# Filter NaN padding out from time and value matrices
time_54_nans = df[('54', 'FRET', 'TIME')].values
time_54 = time_54_nans[~np.isnan(time_54_nans)]

data_54 = np.zeros((1, 2, len(time_54)))
data_54[0, 0, :] = df[('54', 'NBD', 'VALUE')].values[~np.isnan(time_54_nans)]
data_54[0, 1, :] = df[('54', 'FRET', 'VALUE')].values[~np.isnan(time_54_nans)]

time_126 = df[('126', 'FRET', 'TIME')].values
data_126 = np.zeros((1, 2, len(time_126)))
data_126[0, 0, :] = df[('126', 'NBD', 'VALUE')].values
data_126[0, 1, :] = df[('126', 'FRET', 'VALUE')].values

