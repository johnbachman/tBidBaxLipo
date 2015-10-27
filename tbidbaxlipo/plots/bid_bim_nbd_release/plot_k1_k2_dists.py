#from tbidbaxlipo.data.parse_bid_bim_nbd_release import nbd_residues
import cPickle
from itertools import product
import numpy as np
from matplotlib import pyplot as plt
import sys
import re
import os
from tbidbaxlipo.util import format_axis, set_fig_params_for_publication
from tbidbaxlipo.plots.bid_bim_nbd_release.preprocess_data \
        import nbd_residues, max_nbd_value, df

def get_matrix_dimensions(num_nbd_residues, num_reps, num_bin_edges):
    # The number of columns (mutants * reps) * spaces,
    # where spaces = mutants - 1
    num_spaces = num_nbd_residues - 1
    num_cols = (num_nbd_residues * num_reps) + num_spaces
    return (num_bin_edges, num_cols)

def assemble_density_matrix(p_name, filelist):
    # To avoid having to hard-code which residues, activators, and models
    # used, we first iterate over all the provided files and figure out how
    # many things we have in each category (residues, conformations, etc.).
    # To figure out how many unique things we have, we use sets:
    activators = set()
    nbd_residues = set()
    replicates = set()
    pattern = re.compile('pt_data1_(\w+)_NBD_(\d+)_r(\d)_3confs')

    file_dict = {}

    for filename in filelist:
        # First, split off the dirname
        basename = os.path.basename(filename)
        # Next, split off the extension(s)
        prefix = basename.split('.')[0]
        # Next, split the filename into parts at underscores
        m = pattern.match(prefix)
        if not m:
            raise Exception('Could not match filename %s' % prefix)
        # Get the keys from the regex
        (activator, residue, repnum) = m.groups()
        repnum = int(repnum)
        # Load the file
        arr = np.loadtxt(filename)
        # Store the tuple in a dict
        file_dict[(activator, residue, repnum)] = arr
        # Store the keys
        activators.add(activator)
        nbd_residues.add(residue)
        replicates.add(repnum)

    activators = list(sorted(activators))
    nbd_residues = list(sorted(nbd_residues, key=lambda key: int(key)))
    replicates = list(sorted(replicates))

    for act in activators:
        # Sort of a hack so that I don't have to store the setup somewhere else
        if p_name in ['k1', 'k2']:
            num_bin_edges = 51
        elif p_name in ['c1', 'c2']:
            num_bin_edges = 46
        else:
            print "Unknown parameter name: %s" % p_name
            sys.exit(1)

        # Matrix to fill in
        mx_dims = get_matrix_dimensions(len(nbd_residues), len(replicates),
                                        num_bin_edges)
        density_mx = np.zeros(mx_dims)

        # The current column in the matrix
        col_ix = 0
        for res in nbd_residues:
            for rep in replicates:
                # Get the histogram array associated with this rep
                norm_counts = file_dict[(act, res, rep)]
                # Fill in this "stripe" of the matrix
                # We don't fill in the last value b/c it's just there to let us
                # get the axis labels correct
                density_mx[:-1, col_ix] = norm_counts
                col_ix += 1
            # Add a space between the reps
            col_ix += 1
        # Save the matrix
        density_mx_filename = '%s_%s_density_mx.txt' % (act, p_name)
        print("Saving %s" % density_mx_filename)
        np.savetxt(density_mx_filename, density_mx)

def stack_matrices(k1_mx_filename, k2_mx_filename):
    # Load the two matrices
    k1_mx = np.loadtxt(k1_mx_filename)
    k2_mx = np.loadtxt(k2_mx_filename)
    assert k1_mx.shape == k2_mx.shape, \
           "k1/k2 matrices must have same size"

    # Create the combined RGB matrix
    rgb_shape = (k1_mx.shape[1]*2, k1_mx.shape[0], 3)
    #rgb_shape = (k1_mx.shape[1]*2, k1_mx.shape[0])
    #rgb_mx = np.zeros(rgb_shape)
    rgb_mx = np.ones(rgb_shape)
    for row_ix in range(k1_mx.shape[1]):
        rgb_mx[2*row_ix, :, 1] -= k1_mx[:, row_ix]
        rgb_mx[2*row_ix, :, 2] -= k1_mx[:, row_ix]
        rgb_mx[(2*row_ix)+1, :, 0] -= k2_mx[:, row_ix]
        rgb_mx[(2*row_ix)+1, :, 2] -= k2_mx[:, row_ix]
        #rgb_mx[2*row_ix, :] = k1_mx[:, row_ix]
        #rgb_mx[(2*row_ix)+1, :] = k2_mx[:, row_ix]

    return rgb_mx

def merge_matrices(k1_mx_filename, k2_mx_filename, rgb_inverted=True):
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

def plot_matrix(plot_type, density_mx, residues, output_file_base):
    if plot_type == 'k1k2':
        lbound = -6
        ubound = -1
        lbound_ix = 15
        xlabel = 'log$_{10}$(k1, k2)'
    elif plot_type == 'c1c2':
        lbound = 2.5
        ubound = 5.5
        lbound_ix = 0
        xlabel = 'log$_{10}$($F_1$, $F_2$)'
    else:
        raise ValueError("Unknown plot type: %s" % plot_type)

    fig = plt.figure(figsize=(2, 3), dpi=300)
    ax = fig.gca()
    ncols = density_mx.shape[0]
    num_reps = 3
    #assert ncols == (len(residues) * num_reps) + len(residues) - 1, \
    #       "Dimensions of density_mx must match numbers of reps/residues"
    ax.imshow(density_mx[:,lbound_ix:,:], interpolation='none',
              extent=(lbound, ubound, ncols, 0), aspect='auto')
    #ax.matshow(density_mx[:,lbound_ix:], interpolation='none',
    #          extent=(lbound, ubound, ncols, 0), aspect='auto')
    # Plot line representing max observed value
    if plot_type == 'c1c2':
        plt.vlines(np.log10(max_nbd_value), ncols, 0, color='gray')
    # Lines separating the different mutants
    lines = np.arange((2*num_reps), ncols, (2*num_reps) + 1) + 0.5
    #plt.hlines(lines, lbound, ubound)
    """
    # Plot lines representing starting values
    # FIXME FIXME FIXME
    mut_ix_bounds = [0] + list(lines) + [ncols]
    activator = 'Bid'
    row_ix = 0
    for nbd_ix, nbd_residue in enumerate(residues):
        for rep_num in range(1, num_reps + 1):
            timecourse = \
                    df[(activator, 'NBD', nbd_residue, rep_num, 'VALUE')].values
            f0 = timecourse[0]
            #res_ix_lbound = mut_ix_bounds[nbd_ix]
            #res_ix_ubound = mut_ix_bounds[nbd_ix + 1]
            plt.vlines(np.log10(f0), row_ix, row_ix + 1, color='gray')
            row_ix += 1
        row_ix += 1
    """
    # Ticks for the different mutants
    #spacing = 7
    #ytick_positions = np.arange(3, ncols, spacing) + 0.5
    #ax.set_yticks(ytick_positions)
    #assert len(residues) == len(ytick_positions)
    #ax.set_yticklabels(residues)
    # Ticks for the units on the x-axis
    xtick_positions = np.arange(lbound, ubound + 0.1, 0.5)
    ax.set_xticks(xtick_positions)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('NBD Position')
    format_axis(ax)
    fig.subplots_adjust(left=0.17, bottom=0.11, top=0.94)
    fig.savefig('%s.pdf' % output_file_base)
    #fig.savefig('%s.png' % output_file_base)
    import ipdb; ipdb.set_trace()

if __name__ == '__main__':

    usage =  'Usage:\n'
    usage += '    python plot_k1_k2_dists.py assemble [k1|k2|c1|c2] ' \
                                    '[file glob]\n'
    usage += '    python plot_k1_k2_dists.py plot [k1k2|c1c2] ' \
                       'density_file1 density_file2 output_file_base\n'

    if len(sys.argv) < 4:
        print(usage)
        sys.exit(1)
    if sys.argv[1] == 'assemble':
        p_name = sys.argv[2]
        filelist = sys.argv[3:]
        assemble_density_matrix(p_name, filelist)
    elif sys.argv[1] == 'plot' and len(sys.argv) == 6:
        plot_type = sys.argv[2]
        density_file1 = sys.argv[3]
        density_file2 = sys.argv[4]
        output_file_base = sys.argv[5]
        # Check for a correctly specified plot type
        if plot_type not in ['k1k2', 'c1c2']:
            print(usage)
            sys.exit(1)
        # Load and combine the matrices
        #density_mx = merge_matrices(density_file1, density_file2,
        #                              rgb_inverted=True)
        density_mx = stack_matrices(density_file1, density_file2)
        # Get the residues that we're plotting for
        residues = [res for res in nbd_residues if res != 'WT']
        # Plot
        plot_matrix(plot_type, density_mx, residues, output_file_base)
    else:
        print(usage)
        sys.exit(1)

