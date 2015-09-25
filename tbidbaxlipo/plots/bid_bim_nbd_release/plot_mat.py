import pickle
from matplotlib import pyplot as plt
from tbidbaxlipo.util import format_axis, set_fig_params_for_publication
import numpy as np
import sys
from tbidbaxlipo.data.parse_bid_bim_nbd_release import nbd_residues


def combine_matrices(k1_mx_filename, k2_mx_filename, rgb_inverted=True):
    # Load the two matrices
    k1_mx = np.loadtxt(k1_mx_filename)
    k2_mx = np.loadtxt(k2_mx_filename)
    assert k1_mx.shape == k2_mx.shape, \
           "k1/k2 matrices must have same size"

    # Create the combined RGB matrix
    rgb_shape = (k1_mx.shape[1], k1_mx.shape[0], 3)
    if rgb_inverted:
        rgb_mx = np.ones(rgb_shape)
        # k1 goes in the red slice
        rgb_mx[:, :, 1] -= k1_mx.T
        rgb_mx[:, :, 2] -= k1_mx.T
        # k2 goes in the green slice
        rgb_mx[:, :, 0] -= k2_mx.T
        rgb_mx[:, :, 2] -= k2_mx.T
    else:
        rgb_mx = np.zeros(rgb_shape)
        # k1 goes in the red slice
        rgb_mx[:, :, 0] = k1_mx.T
        # k1 goes in the green slice
        rgb_mx[:, :, 1] = k2_mx.T
    return rgb_mx

if __name__ == '__main__':
    set_fig_params_for_publication()

    residues = [res for res in nbd_residues if res != 'WT']

    #k1_filename = 'mcmc/Bid_k1_density_mx.txt'
    #k2_filename = 'mcmc/Bid_k2_density_mx.txt'
    #k1_filename = 'mcmc/Bid_c1_density_mx.txt'
    #k2_filename = 'mcmc/Bid_c2_density_mx.txt'

    plt.ion()
    density_mx = combine_matrices(k1_filename, k2_filename,
                                  rgb_inverted=True)

    #density_filename = sys.argv[1]
    #density_mx = np.loadtxt(density_filename)
    #inverted_mx = np.ones(density_mx.shape)
    #inverted_mx = inverted_mx - density_mx

    ncols = density_mx.shape[0]
    #lbound = -6
    #ubound = -1
    lbound = 2.5
    ubound = 5.5
    fig = plt.figure(figsize=(2, 3), dpi=300)
    ax = fig.gca()
    #ax.matshow(inverted_mx.T, cmap='gray')
    #, extent=(-4, -1, 0, 19), aspect='auto')
    ax.imshow(density_mx[:,15:,:], interpolation='none',
              extent=(lbound, ubound, ncols, 0), aspect='auto')
    format_axis(ax)
    lines = np.arange(3, ncols, 4) + 0.5

    plt.hlines(lines, lbound, ubound)

    ytick_positions = np.arange(1, ncols, 4) + 0.5
    ax.set_yticks(ytick_positions)
    assert len(residues) == len(ytick_positions)
    ax.set_yticklabels(residues)
    #xtick_positions = np.arange(-4, -0.9, 0.5)
    #ax.set_xticks(xtick_positions)
    ax.set_xlabel('log$_{10}$(k)')
    ax.set_ylabel('NBD Position')
    fig.subplots_adjust(left=0.17, bottom=0.11, top=0.94)
    fig.savefig('k1_k2_posteriors.pdf')
    #plt.show()

