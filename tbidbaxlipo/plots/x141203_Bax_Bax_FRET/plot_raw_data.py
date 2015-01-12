import numpy as np
from matplotlib import pyplot as plt
from remove_outliers import df
from parse_data import nbd_residues

dtype_line_colors = {'Release':'r', 'NBD':'g', 'FRET':'b'}

def plot_all(df):
    for mutant in nbd_residues:
        plt.figure()

        for dtype in ['Release', 'FRET', 'NBD']:
            #t = df[(dtype, mutant, 'TIME')]
            v = df[(dtype, mutant, 'VALUE')]
            norm_v = v / np.max(v)
            #plt.plot(t, v, dtype_line_colors[dtype])
            plt.plot(norm_v, dtype_line_colors[dtype])
        plt.xlabel('Time (sec)')
        plt.ylabel('Norm Signal')
        plt.title(mutant)

if __name__ == '__main__':
    plt.ion()
    plot_all(df)
