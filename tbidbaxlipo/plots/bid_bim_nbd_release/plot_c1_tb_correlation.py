# Load the data from the CSV files

import csv
from matplotlib import pyplot as plt
from itertools import product
from scipy.stats import linregress
import numpy as np

rel_dict = {}
c1_dict = {}
activators = set()
nbd_residues = set()
replicates = set()

with open('data1_release_peak_times.csv') as rel_file:
    rel_reader = csv.reader(rel_file, delimiter=',')
    for row in rel_reader:
        (act, res, rep, time) = row
        rel_dict[(act, res, rep)] = time
        activators.add(act)
        nbd_residues.add(res)
        replicates.add(rep)

with open('mcmc/pt_data1_c1_peak_times.csv') as c1_file:
    c1_reader = csv.reader(c1_file, delimiter=',')
    for row in c1_reader:
        (act, res, rep, avg, sd) = row
        c1_dict[(act, res, rep)] = (avg, sd)
        activators.add(act)
        nbd_residues.add(res)
        replicates.add(rep)

activators = list(sorted(activators))
nbd_residues = list(sorted(nbd_residues))
replicates = list(sorted(replicates))

# Plot
plt.ion()
xdata = []
ydata = []

ctrl = np.linspace(0, 1000, 5)
plt.figure()
plt.plot(ctrl, ctrl, color='gray')
for (act, nbd, rep) in product(activators, nbd_residues, replicates):
    if nbd == 'WT': # or nbd == '68' or nbd == '79' or nbd == '120' or \
        continue
    #   nbd =='188':
    cur = None
    #if nbd == '54' or nbd == '62':
    #    pass
    #else:
    #    continue
    rel_time = rel_dict[(act, nbd, rep)]
    (c1_avg, c1_sd) = c1_dict[(act, nbd, rep)]
    (rel_time, c1_avg, c1_sd) = map(float, [rel_time, c1_avg, c1_sd])
    xdata.append(rel_time)
    ydata.append(c1_avg)
    #print act, nbd, rep
    #print rel_time, float(c1_avg), float(c1_sd)
    #plt.errorbar(float(rel_time), float(c1_avg), yerr=float(c1_sd))
    if nbd == cur:
        col = 'r'
    else:
        col = 'b'
    plt.scatter(rel_time, c1_avg, color=col)
    print act, nbd, rep
    #import ipdb; ipdb.set_trace()

xdata = np.array(xdata)
ydata = np.array(ydata)
(slope, intercept, r, p, stderr) = linregress(xdata, ydata)

plt.plot(ydata, ydata*slope + intercept)
print slope, intercept, r, p

