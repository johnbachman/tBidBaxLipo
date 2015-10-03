from tbidbaxlipo.data.parse_bid_bim_fret_nbd_release import df, nbd_residues
from preprocess_data import df_pre
from itertools import product
import numpy as np
from matplotlib import pyplot as plt

plt.ion()

"""
for rep in [1, 2, 3]:
    dkey = ('Bim', 'FRET', '184', rep, 'VALUE')
    v = df[dkey].values
    vpre = df_pre[dkey].values
    ind = np.arange(len(v))
    plt.subplot(3, 1,rep)
    plt.plot(v, linestyle='', marker='.', color='b', markersize=10)
    plt.plot(vpre, color='red', marker='.', linestyle='', markersize=10)
    plt.hlines(0, 0, len(v))
"""

def plot_fret_outliers(df, df_pre, activators, residues, reps=(1, 2, 3)):
    for act in activators:
        for res in residues:
            plt.figure('%s, NBD-%sC-Bax' % (act, res), figsize=(6, 8))
            for rep in reps:
                plt.subplot(3, 1, rep)
                plt.title('%s, NBD-%sC-Bax, rep %d' % (act, res, rep))
                key = (act, 'FRET', res, rep, 'VALUE')
                v = df[key].values
                vpre = df_pre[key].values
                plt.plot(v, linestyle='', marker='.', color='b', markersize=10)
                plt.plot(vpre, color='red', marker='.', linestyle='',
                         markersize=10)
                plt.ylabel('% FRET')
            # Add a line at 0
            plt.hlines(0, 0, len(v))
            # Add axis labels
            plt.xlabel('Timepoint index')
            plt.subplots_adjust(hspace=0.4, top=0.93)

if __name__ == '__main__':
    plot_fret_outliers(df, df_pre, ('Bid',), ('62', '126'))

