"""A script that iterates over a series of .mcmc files from fits of the
three-conformation model, loads the chains and calculates the average value
(and SD) associated with the peak of C1 intermediate, at
log(k1/k2)/(k1 - k2)."""

import sys
import os.path
import pickle
import numpy as np

filename = sys.argv[1]

print("Loading: %s" % filename)
basename = os.path.basename(filename)
with open(filename) as f:
    (gf, sampler) = pickle.load(f)

# The index of k1 in the chain is the same as the index in estimate_params
k1 = gf.builder.model.parameters['c0_to_c1_k']
k2 = gf.builder.model.parameters['c1_to_c2_k']
k1_index = gf.builder.estimate_params.index(k1)
k2_index = gf.builder.estimate_params.index(k2)
k1_samples = 10 ** sampler.flatchain[0, :, k1_index]
k2_samples = 10 ** sampler.flatchain[0, :, k2_index]
times = (np.log(k1_samples / k2_samples) /
         (k1_samples - k2_samples))
time_avg = np.mean(times)
time_sd = np.std(times, ddof=1)
print time_avg
print time_sd

#from matplotlib import pyplot as plt
#plt.ion()
#plt.hist(times, bins=20)

output_filename = '%s.c1_peak' % basename
print("Writing: %s" % output_filename)
with open(output_filename, 'w') as f:
    f.write('%s, %s' % (time_avg, time_sd))

