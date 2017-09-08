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

def get_parameter_samples(pname, gf, sampler):
    # The index of p in the chain is the same as the index in estimate_params
    p = gf.builder.model.parameters[pname]
    p_index = gf.builder.estimate_params.index(p)
    p_samples = sampler.flatchain[0, :, p_index]
    return p_samples

def parameter_histogram(p_samples, lower_bound, upper_bound, filename, p_scale,
                        num_bins):
    # -- Calculate the histograms for the parameter (normalized density) --
    # For now, we declare the width of the prior up front to avoid having to
    # load it first and then set it
    width = upper_bound - lower_bound
    # Get bin boundaries for prior distribution
    #num_bin_edges = width * bins_per_log + 1
    num_bin_edges = num_bins + 1
    bin_edges = np.linspace(lower_bound, upper_bound, num_bin_edges)
    # Calculate the histogram
    if p_scale == 'lin':
        p_samples = 10 ** p_samples
    (p_counts, _) = np.histogram(p_samples, bins=bin_edges)
    p_norm_counts = p_counts / float(len(p_samples))
    # Save histograms
    print("Writing: %s" % filename)
    np.savetxt(filename, p_norm_counts)

def c1_peak_times(k1_log_samples, k2_log_samples, basename):
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

    # Histograms
    param_info = [('c0_to_c1_k', 'k1_hist', -6, -1, 50, 'log'),
                  ('c1_to_c2_k', 'k2_hist', -6, -1, 50, 'log'),
                  ('c1_scaling', 'c1_scaling', 0, 10, 50, 'lin', ),
                  ('c2_scaling', 'c2_scaling', 0, 10, 50, 'lin'),
                  ('fret1_scaling', 'fret1_scaling', 0, 100, 50, 'lin'),
                  ('fret2_scaling', 'fret2_scaling', 0, 100, 50, 'lin')]
    for (p_name, file_suffix, lower_bound, upper_bound, bins, p_scale) \
                                                            in param_info:
        p_samples = get_parameter_samples(p_name, gf, sampler)
        parameter_histogram(p_samples, lower_bound, upper_bound,
                            '%s.%s' % (filename, file_suffix), p_scale,
                            num_bins=bins)
    # Peak times
    k1_samples = get_parameter_samples('c0_to_c1_k', gf, sampler)
    k2_samples = get_parameter_samples('c1_to_c2_k', gf, sampler)
    c1_peak_times(k1_samples, k2_samples, filename)

