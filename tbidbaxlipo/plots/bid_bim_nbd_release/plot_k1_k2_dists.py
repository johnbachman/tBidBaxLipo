#from tbidbaxlipo.data.parse_bid_bim_nbd_release import nbd_residues
import cPickle
from itertools import product
import numpy as np
from matplotlib import pyplot as plt

activator = 'Bid'
nbd_residues = ['3', '5', '15', '54', '68']
reps = [1, 2, 3]
confs = 3

# Create a matrix with bins of the distribution in the y axis and entries
# for the mutants/reps on the x-axis.

# For now, we declare the width of the prior up front to avoid having to
# load it first and then set it
upper_bound = -1
lower_bound = -4
width = upper_bound - lower_bound
# Get bin boundaries for prior distribution
bins_per_log = 10
num_bin_edges = width * bins_per_log + 1
bin_edges = np.linspace(lower_bound, upper_bound, num_bin_edges)
num_bins = num_bin_edges - 1
print num_bins
# The number of columns (mutants * reps) * spaces,
# where spaces = mutants - 1
num_spaces = len(nbd_residues) - 1
num_cols = (len(nbd_residues) * len(reps)) + num_spaces
# Matrix to fill in
density_mx = np.ones((num_bin_edges, num_cols))

# The current column in the matrix
col_ix = 0
for res in nbd_residues:
    for rep in reps:
        filename = 'mcmc/pt_data1_%s_NBD_%s_r%s_%sconfs.mcmc' % \
                   (activator, res, rep, confs)
        # Load the data
        print("Loading %s" % filename)
        with open(filename) as f:
            (gf, sampler) = cPickle.load(f)
        # Get the index for parameter k1
        p = gf.builder.model.parameters['c0_to_c1_k']
        p_ix = gf.builder.estimate_params.index(p)
        # Get the samples for this parameter
        p_samples = sampler.flatchain[0, :, p_ix]
        num_samples = len(p_samples)
        print np.mean(p_samples)
        # Get bounds for the prior
        #p_prior = gf.builder.priors[p_ix]
        (counts, _) = np.histogram(p_samples, bins=bin_edges)
        norm_counts = counts / float(num_samples)
        # Fill in this "stripe" of the matrix
        # We don't fill in the last value b/c it's just there to let us
        # get the axis labels correct
        density_mx[:-1, col_ix] = 1 - norm_counts
        col_ix += 1
    # Add a space between the reps
    col_ix += 1

print density_mx

# Get a colormap
#plt.get_cmap('gray')
#plt.imshw
plt.ion()
