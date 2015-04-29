"""A script that iterates over a series of .mcmc files, loads the chains, gets
the maximum a posteriori probability and prints the sorted list."""

import sys
import os.path
import pickle
import numpy as np

file_list = sys.argv[1:]
maxp_list = []
evidence_list = []
for filename in file_list:
    print("Loading: %s" % filename)
    basename = os.path.basename(filename)
    with open(filename) as f:
        (gf, sampler) = pickle.load(f)

    lnprob = sampler.lnprobability[0]
    maxp_ix = np.unravel_index(np.argmax(lnprob), lnprob.shape)
    maxp = lnprob[maxp_ix[0], maxp_ix[1]]
    maxp_list.append((basename, maxp))
    evidence = sampler.thermodynamic_integration_log_evidence()
    evidence_list.append((basename, evidence))

# Sort on the probability
sorted_maxp = sorted(maxp_list, key=lambda tup: tup[1], reverse=True)
sorted_evidence = sorted(evidence_list, key=lambda tup: tup[1], reverse=True)
print("--- Max a posteriori: ---")
for (filename, maxp) in sorted_maxp:
    print("%s: %s" % (maxp, filename))
print("\n")
print("--- log(evidence) ---")
for (filename, evidence) in sorted_evidence:
    print("%s: %s" % (evidence, filename))
