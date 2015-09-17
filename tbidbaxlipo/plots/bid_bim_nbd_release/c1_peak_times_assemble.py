"""Iterate all of the filenames containing the extracted c1 peak
values from MCMC runs (provided as a glob at the command line) and tabulate
them as a single CSV file and pandas dataframe.

Example usage::

    python -m tbidbaxlipo.plots.bid_bim_nbd_release.c1_peak_times_assemble *_3confs.mcmc.c1_peak
"""

import sys
import re
import csv
from itertools import product
import numpy as np

filelist = sys.argv[1:]

# To avoid having to hard-code which residues, activators, and models used,
# we first iterate over all the provided files and figure out how many things
# we have in each category (residues, conformations, etc.).
# To figure out how many unique things we have, we use sets:
activators = set()
nbd_residues = set()
replicates = set()
pattern = re.compile('pt_data1_(\w+)_NBD_(\d+)_r(\d)_3confs')

c1_dict = {}

for filename in filelist:
    # First, split off the extension(s)
    prefix = filename.split('.')[0]
    # Next, split the filename into parts at underscores
    m = pattern.match(prefix)
    if not m:
        raise Exception('Could not match filename %s' % prefix)
    # Get the keys from the regex
    (activator, residue, repnum) = m.groups()
    repnum = int(repnum)
    # Load the file
    arr = np.loadtxt(filename, delimiter=',')
    # Store the tuple in a dict
    c1_dict[(activator, residue, repnum)] = arr
    # Store the keys
    activators.add(activator)
    nbd_residues.add(residue)
    replicates.add(repnum)

activators = list(sorted(activators))
nbd_residues = list(sorted(nbd_residues, key=lambda key: int(key)))
replicates = list(sorted(replicates))

# CSV Writer
with open('pt_data1_c1_peak_times.csv', 'w') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    for (act, nbd, rep) in product(activators, nbd_residues, replicates):
        (time_avg, time_sd) = c1_dict[(act, nbd, rep)]
        csvwriter.writerow([act, nbd, rep, time_avg, time_sd])

