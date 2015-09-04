"""A script that iterates over a series of .mcmc files, loads the chains, gets
the maximum a posteriori probability and prints the sorted list."""

import sys
import os.path
import pickle
#import numpy as np

filename = sys.argv[1]

print("Loading: %s" % filename)
basename = os.path.basename(filename)
with open(filename) as f:
    (gf, sampler) = pickle.load(f)

evidence = sampler.thermodynamic_integration_log_evidence()
output_filename = '%s.evi' % basename
print("Writing: %s" % output_filename)
with open(output_filename, 'w') as f:
    f.write('%s, %s' % evidence)
