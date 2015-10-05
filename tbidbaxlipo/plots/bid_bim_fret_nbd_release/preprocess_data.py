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

# Bid/47/1: two early outliers
outliers = [7, 15]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '47', 1, 'VALUE')][outlier_ix] = np.nan
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

# Bid/122/1
outliers = [9]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '122', 1, 'VALUE')][outlier_ix] = np.nan
# Bid/122/2: pretty good
# Bid/122/3: pretty good

# Bim/122/1: pretty good
# Bim/122/2
outliers = [7, 8, 28, 45, 84, 89, 90, 95, 96]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '122', 2, 'VALUE')][outlier_ix] = np.nan
# Bim/122/3
outliers = [70, 71, 84, 86, 96]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '122', 3, 'VALUE')][outlier_ix] = np.nan

# 126
# Bid/126/1: OK
# Bid/126/2: OK
# Bid/126/3: OK

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

# Bid/138/1
outliers = [7, 13, 14, 36, 47, 64, 81, 99, 100, 106, 107, 110, 113]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '138', 1, 'VALUE')][outlier_ix] = np.nan
# Bid/138/2: this one in rough shape, many drops in signal, hard
# to draw a line for what to take out
outliers = [2, 19, 24, 40, 43, 47, 52, 74, 80, 81, 82, 85, 86, 89, ]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '138', 2, 'VALUE')][outlier_ix] = np.nan
# Bid/138/3
outliers = [7, 8, 9, 13, 31, 42, 51, 54, 55, 58, 73]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '138', 3, 'VALUE')][outlier_ix] = np.nan

# Bim/138/1
outliers = [53, 56, 62, 75, 90, 103, 114]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '138', 1, 'VALUE')][outlier_ix] = np.nan
# Bim/138/2
outliers = [41]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '138', 2, 'VALUE')][outlier_ix] = np.nan
# Bim/138/3
outliers = [32]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '138', 3, 'VALUE')][outlier_ix] = np.nan

# Bid/151/1
outliers = [2, 6, 76, 83, 88, 99, 103, 106, 118]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '151', 1, 'VALUE')][outlier_ix] = np.nan
# Bid/151/2: noise, but no outliers
# Bid/151/3: noise, but no outliers

# Bim/151/1
outliers = [10]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '151', 1, 'VALUE')][outlier_ix] = np.nan
# Bim/151/2
outliers = [1, 50, 64, 87, 89, 94, 101, 103]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '151', 2, 'VALUE')][outlier_ix] = np.nan
# Bim/151/3: noisy, but no major outliers

# Bid/175/1: very noisy, but hard to pick out outliers
outliers = [32, 41, 42, 51, 93, 107]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '175', 1, 'VALUE')][outlier_ix] = np.nan
# Bid/175/2: pretty good
# Bid/175/3: pretty good

# Bim/175/1: not too bad
# Bim/175/2
outliers = [55, 57, 68, 71, 75, 87, 88, 106, 108, 114]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '175', 2, 'VALUE')][outlier_ix] = np.nan
# Bim/175/3
outliers = [44, 56, 57, 67, 69, 70, 72, 78]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '175', 3, 'VALUE')][outlier_ix] = np.nan

# Bid/184/1: pretty good
# Bid/184/2
outliers = [26, 49, 52, 85, 94]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '184', 2, 'VALUE')][outlier_ix] = np.nan
# Bid/184/3
outliers = [2, 4, 39]
for outlier_ix in outliers:
    df_pre[('Bid', 'FRET', '184', 3, 'VALUE')][outlier_ix] = np.nan

# Bim/184/1: pretty good
# Bim/184/2: pretty good
# Bim/184/3: pretty good


