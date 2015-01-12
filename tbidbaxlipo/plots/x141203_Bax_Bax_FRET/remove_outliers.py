from preprocess_data import df
import numpy as np

# There are a number of outliers in the FRET data that Justin attributes
# to debris floating in the cuvette that produces momentary spikes. Here
# we explicitly remove several of these from the 54C data, setting them to
# NaN.
outliers = [21, 24, 77, 103, 108, 113]
for outlier_ix in outliers:
    df[('FRET', '54', 'VALUE')][outlier_ix] = np.nan



