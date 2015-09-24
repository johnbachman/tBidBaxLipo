"""A script that iterates over a series of .mcmc files from fits of the
three-conformation model, loads the chains and calculates the average value
(and SD) associated with the peak of C1 intermediate, at
log(k1/k2)/(k1 - k2).
"""

import sys
import os.path
import pickle
import numpy as np

def load_mcmc_file(filename):
    print("Loading: %s" % filename)
    basename = os.path.basename(filename)
    with open(filename) as f:
        (gf, sampler) = pickle.load(f)
    return (gf, sampler)

def get_k1_k2_samples(gf, sampler):
    # The index of k1 in the chain is the same as the index in estimate_params
    k1 = gf.builder.model.parameters['c0_to_c1_k']
    k2 = gf.builder.model.parameters['c1_to_c2_k']
    k1_index = gf.builder.estimate_params.index(k1)
    k2_index = gf.builder.estimate_params.index(k2)
    k1_log_samples = sampler.flatchain[0, :, k1_index]
    k2_log_samples = sampler.flatchain[0, :, k2_index]
    return (k1_log_samples, k2_log_samples)

def k1_k2_histogram(k1_samples, k2_samples):
    # -- Calculate the histograms for k1 and k2 (normalized) --
    # For now, we declare the width of the prior up front to avoid having to
    # load it first and then set it
    upper_bound = -1
    lower_bound = -6
    width = upper_bound - lower_bound
    # Get bin boundaries for prior distribution
    bins_per_log = 10
    num_bin_edges = width * bins_per_log + 1
    bin_edges = np.linspace(lower_bound, upper_bound, num_bin_edges)
    num_bins = num_bin_edges - 1
    # Calculate the histogram
    (k1_counts, _) = np.histogram(k1_samples, bins=bin_edges)
    (k2_counts, _) = np.histogram(k2_samples, bins=bin_edges)
    k1_norm_counts = k1_counts / float(len(k1_samples))
    k2_norm_counts = k2_counts / float(len(k2_samples))
    # Save histograms
    k1_hist_filename = '%s.k1_hist' % basename
    k2_hist_filename = '%s.k2_hist' % basename
    print("Writing %s" % k1_hist_filename)
    print("Writing %s" % k2_hist_filename)
    np.savetxt(k1_hist_filename, k1_norm_counts)
    np.savetxt(k2_hist_filename, k2_norm_counts)

def c1_peak_times(k1_log_samples, k2_log_samples):
    # Transform the samples to linear space
    k1_samples = 10 ** k1_log_samples
    k2_samples = 10 ** k2_log_samples

    # -- Calculate and save the peak times --
    times = (np.log(k1_samples / k2_samples) /
             (k1_samples - k2_samples))
    time_avg = np.mean(times)
    time_sd = np.std(times, ddof=1)

    peak_time_filename = '%s.c1_peak' % basename
    print("Writing: %s" % peak_time_filename)
    with open(peak_time_filename, 'w') as f:
        f.write('%s, %s' % (time_avg, time_sd))

if __name__ == '__main__':
    filename = sys.argv[1]

    (gf, sampler) = load_mcmc_file(filename)
    # if three conformations:
    (k1_samples, k2_samples) = get_k1_k2_samples(gf, sampler)
    k1_k2_histogram(k1_samples, k2_samples)
    c1_peak_times(k1_samples, k2_samples)

