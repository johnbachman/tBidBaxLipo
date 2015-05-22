import numpy as np
from matplotlib import pyplot as plt
import sys

if len(sys.argv) < 2:
    print("Specify the plot filename.")
    sys.exit()

plot_filename = sys.argv[1]

num_entries = 15
num_reps = 3
results = np.zeros((num_entries, num_reps))
num_procs_list = range(2, 17)
for ix, num_procs in enumerate(num_procs_list):
    filename = '%s.txt' % num_procs
    with open(filename) as f:
        lines = f.readlines()
        for rep_ix, line in enumerate(lines):
            results[ix, rep_ix] = line.strip()

m = np.mean(results, axis=1)
sd = np.std(results, axis=1, ddof=1)
plt.ion()
plt.errorbar(num_procs_list, m, yerr=sd)
plt.xlabel("Num processes")
plt.ylabel("Time, 100 MCMC steps (sec)")
plt.savefig('%s.pdf' % plot_filename)
