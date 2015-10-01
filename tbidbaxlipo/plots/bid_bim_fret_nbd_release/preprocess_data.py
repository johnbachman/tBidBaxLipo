import sys
import numpy as np
from tbidbaxlipo.data.parse_bid_bim_fret_nbd_release import df, nbd_residues
from copy import deepcopy

# There are a number of outliers in the FRET data that Justin attributes
# to debris floating in the cuvette that produces momentary spikes. Here
# we explicitly remove several of these from the 54C data, setting them to
# NaN.
df_pre = deepcopy(df)

# Bid/3/1: though negative, no terrible outliers
# Bid/3/2: no terrible outliers
# Bid 3/3: no terrible outliers

# Bim/3/1
outliers = [55, 58, 94, 105]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '3', 1, 'VALUE')][outlier_ix] = np.nan
# Bim/3/2
outliers = [12, 18, 20, 21, 30, 48, 56, 61, 69, 100, 106, 111, 114]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '3', 2, 'VALUE')][outlier_ix] = np.nan
# Bim/3/3: not too bad

# Bid/15/1: goes negative, but not bad
# Bid/15/2: not too bad
# Bid/15/3: not too bad

# Bim/15/1
outliers = [39, 41, 43, 61, 76, 84]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '15', 1, 'VALUE')][outlier_ix] = np.nan
# Bim/15/2
outliers = [3, 7, 8, 11, 13]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '15', 2, 'VALUE')][outlier_ix] = np.nan
# Bim/15/3
outliers = [84]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '15', 3, 'VALUE')][outlier_ix] = np.nan

# Bid/36/1: pretty good
# Bid/36/2: pretty good
# Bid/36/3: pretty good

# Bim/36/1
outliers = [21, 44, 59, 61, 63, 75, 82, 95, 96]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '36', 1, 'VALUE')][outlier_ix] = np.nan
# Bim/36/2
# Bim/36/3
outliers = [0, 28, 68, 97, 100]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '36', 3, 'VALUE')][outlier_ix] = np.nan

# Bid/47/1: OK, ignoring two early outliers though
# Bid/47/2: Several outlying negative pts:
outliers = [86, 97, 110, 111, 116]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '47', 2, 'VALUE')][outlier_ix] = np.nan
# Bid/47/3: No outliers

# Bim/47/1: no outliers
# Bim/47/2
outliers = [1, 13, 20, 22, 24, 25, 54, 55, 57, 65, 68, 70, 71, 72, 76, 97, 117]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '47', 2, 'VALUE')][outlier_ix] = np.nan
# Bim/47/3
outliers = [33, 39, 84]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '47', 3, 'VALUE')][outlier_ix] = np.nan

# Bid/54/1
# Bid/54/2
outliers = [111]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '54', 2, 'VALUE')][outlier_ix] = np.nan
# Bid/54/3

# Bim/54/1: no outliers
# Bim/54/2: no outliers
# Bim/54/3: two bad negative ones, plus several more dubious ones
outliers = [10, 21]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '54', 3, 'VALUE')][outlier_ix] = np.nan

# Bid/62/1: no outliers
# Bid/62/2: many outliers, some worse than others
outliers = [23, 32, 44, 53, 56, 68, 72, 81, 83, 98, 112, 113]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '62', 2, 'VALUE')][outlier_ix] = np.nan
# Bid/62/3: this one is in rough shape. Nearly the entire trajectory is
# negative, with only major outlying spikes above 0. This removes the spikes
# but doesn't deal with the negative problem.
outliers = [0, 1, 3, 11, 12, 15, 20, 29, 31, 33, 36, 37, 41, 44, 51, 53, 61,
            64, 85, 86, 94, 96, 97, 98, 100, 111]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '62', 3, 'VALUE')][outlier_ix] = np.nan

# Bim/62/1
outliers = [110]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '62', 1, 'VALUE')][outlier_ix] = np.nan
# Bim/62/2: this one is almost unusable, with negative values as low as
# -169% FRET! Since most of the consistent points are positive, the simplest
# heuristic is to just throw away the negative values.
outliers = [1, 5, 7, 13, 17, 24, 25, 26, 27, 29, 35, 36, 38,
            42, 43, 46, 47, 48, 50, 51, 53, 55, 60, 65, 66, 69,
            70, 77, 78, 82, 88, 90, 94, 95, 97, 100, 104, 105, 106,
            107, 108, 109, 113, 115, 118]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '62', 2, 'VALUE')][outlier_ix] = np.nan
# Bim/62/3
outliers = [30, 62, 85, 91, 92, 101]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '62', 3, 'VALUE')][outlier_ix] = np.nan

# 122

# 126
# Bim/126/1
outliers = [49]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '126', 1, 'VALUE')][outlier_ix] = np.nan
# Bim/126/2: mostly negative, with few positive outliers
outliers = [ 2, 3, 7, 13, 15, 25, 29, 31, 37, 40, 43, 47, 57, 71, 74,
            76, 78, 82, 102]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '126', 2, 'VALUE')][outlier_ix] = np.nan
# Bim/126/3
outliers = [48, 49]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '126', 3, 'VALUE')][outlier_ix] = np.nan

# 138
# 151
# 175
# 184
