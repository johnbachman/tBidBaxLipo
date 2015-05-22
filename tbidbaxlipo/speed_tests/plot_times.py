import numpy as np
from matplotlib import pyplot as plt
import sys

if len(sys.argv) < 3:
    print("Specify the plot filename and the number of reps.")
    sys.exit()

plot_filename = sys.argv[1]
num_reps = int(sys.argv[2])

num_procs_list = range(16, 65, 4)
#num_procs_list = range(2, 2 + num_entries)
num_entries = len(num_procs_list)
results = np.zeros((num_entries, num_reps))
nworkers_list = [procs - 1 for procs in num_procs_list]
for ix, num_procs in enumerate(num_procs_list):
    filename = '%s.txt' % num_procs
    with open(filename) as f:
        lines = f.readlines()
        for rep_ix, line in enumerate(lines):
            time = float(line.strip())
            results[ix, rep_ix] = time * nworkers_list[ix]

m = np.mean(results, axis=1)
plt.ion()

if num_reps > 1:
    sd = np.std(results, axis=1, ddof=1)
    plt.errorbar(nworkers_list, m, yerr=sd, marker='o')
else:
    plt.plot(nworkers_list, m, marker='o')
    plt.ylim([0, None])
plt.xlabel("Num workers")
plt.ylabel("Normalized 1-worker time, MCMC steps (sec)")
plt.savefig('%s.pdf' % plot_filename)
