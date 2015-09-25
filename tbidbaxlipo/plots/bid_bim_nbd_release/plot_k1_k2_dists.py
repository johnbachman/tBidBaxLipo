#from tbidbaxlipo.data.parse_bid_bim_nbd_release import nbd_residues
import cPickle
from itertools import product
import numpy as np
from matplotlib import pyplot as plt
import sys
import re
import os

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

if __name__ == '__main__':

    usage =  'Usage:\n'
    usage += '    python plot_k1_k2_dists.py assemble [k1|k2] ' \
                                    '[*.k1_hist|*.k2_hist]\n'
    usage += '    python plot_k1_k2_dists.py plot k1_density_file ' \
                                    'k2_density_file\n'
    if len(sys.argv) < 4:
        print(usage)
        sys.exit(1)
    if sys.argv[1] == 'assemble':
        p_name = sys.argv[2]
        filelist = sys.argv[3:]
        assemble_density_matrix(p_name, filelist)
    elif sys.argv[1] == 'plot':
        pass
    else:
        print(usage)
        sys.exit(1)

#activator = 'Bid'
#nbd_residues = ['3', '5', '15', '54', '68']
#reps = [1, 2, 3]
#confs = 3

# Create a matrix with bins of the distribution in the y axis and entries
# for the mutants/reps on the x-axis.

# Get a colormap
#plt.get_cmap('gray')
#plt.imshw
#plt.ion()


