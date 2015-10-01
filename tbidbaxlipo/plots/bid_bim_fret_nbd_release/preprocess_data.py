import sys
import numpy as np
from tbidbaxlipo.data.parse_bid_bim_fret_nbd_release import df, nbd_residues
from copy import deepcopy

# There are a number of outliers in the FRET data that Justin attributes
# to debris floating in the cuvette that produces momentary spikes. Here
# we explicitly remove several of these from the 54C data, setting them to
# NaN.
df_pre = deepcopy(df)

# Outliers in Bim 47C FRET data
outliers = [1, 13, 20, 22, 24, 25, 54, 55, 57, 65, 68, 70, 71, 72, 76, 97, 117]
for outlier_ix in outliers:
    df_pre[('Bim', 'FRET', '47', 2, 'VALUE')][outlier_ix] = np.nan

# 3
"""
df_pre[('Bid', 'FRET', '47', 1, 'VALUE')][outlier_ix] = np.nan
df_pre[('Bim', 'FRET', '47', 1, 'VALUE')][outlier_ix] = np.nan
df_pre[('Bid', 'FRET', '47', 2, 'VALUE')][outlier_ix] = np.nan
df_pre[('Bim', 'FRET', '47', 2, 'VALUE')][outlier_ix] = np.nan
df_pre[('Bid', 'FRET', '47', 3, 'VALUE')][outlier_ix] = np.nan
df_pre[('Bim', 'FRET', '47', 3, 'VALUE')][outlier_ix] = np.nan
"""

# 15
# 36
# 47
# 54
# 62
# 122
# 126
# 138
# 151
# 175
# 184
