from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

font = {'size': 8}
matplotlib.rc('font', **font)

line_colors = {'Bid': 'r', 'Bim': 'b'}
line_styles = {1:':', 2:'-', 3:'--'}

def plot_all():
    for mutant in nbd_residues:
        plt.figure(figsize=(11, 5))
        # Make the release plot
        plt.subplot(1, 2, 1)
        for activator in ['Bid', 'Bim']:
            for i in range(1, 4):
                t = df[(activator, 'Release', mutant, i, 'TIME')]
                v = df[(activator, 'Release', mutant, i, 'VALUE')]
                plt.plot(t, v, label='%s Rep %d' % (activator, i),
                        color=line_colors[activator],
                        linestyle=line_styles[i])

                plt.xlabel('Time (sec)')
                plt.ylabel('Pct. Release')
                plt.title('Release for NBD-%s-Bax' % mutant)
                plt.legend(loc='lower right')
        # There is no NBD curve for WT Bax, so skip the NBD
        # plot
        if mutant == 'WT':
            continue
        # Make the NBD plot
        plt.subplot(1, 2, 2)
        for activator in ['Bid', 'Bim']:
            for i in range(1, 4):
                t = df[(activator, 'NBD', mutant, i, 'TIME')]
                v = df[(activator, 'NBD', mutant, i, 'VALUE')]
                plt.plot(t, v, label='%s Rep %d' % (activator, i),
                        color=line_colors[activator],
                        linestyle=line_styles[i])
                plt.xlabel('Time (sec)')
                plt.ylabel('F/F0')
                plt.title('F/F0 for NBD-%s-Bax' % mutant)
                plt.legend(loc='lower right')
        plt.tight_layout()
        plt.show()

rt = df[('Bid', 'Release', '68', 3, 'TIME')].values
ry = df[('Bid', 'Release', '68', 3, 'VALUE')].values
nt = df[('Bid', 'NBD', '68', 3, 'TIME')].values
ny = df[('Bid', 'NBD', '68', 3, 'VALUE')].values

import titration_fits as tf

plt.ion()

plt.figure()
plt.plot(rt, ry)

twoexp = tf.TwoExp()
params = twoexp.fit_timecourse(rt, ry)

plt.plot(rt, twoexp.fit_func(rt, params))

# Pores
#plt.figure()
#pores = -np.log(1 - (ry/100.))
#plt.plot(rt, pores)

# NBD
#twoexp = tf.TwoExp()
#params = twoexp.fit_timecourse(nt, ny)
#plt.plot(nt, twoexp.fit_func(nt, params))

from tbidbaxlipo.models.nbd.multiconf import Builder
from pysb.integrate import Solver
from tbidbaxlipo.util import fitting
b = Builder()
b.build_model_multiconf(3, ny[0], normalized_data=True)
Bax_0 = 1.
c0_scaling = fitting.Parameter(1)
c0_to_c1_k = fitting.Parameter(2e-3)
c1_scaling = fitting.Parameter(2.)
c1_to_c2_k = fitting.Parameter(1e-4)
c2_scaling = fitting.Parameter(1)

s = Solver(b.model, nt)
def model_func(t):
    s.run(param_values=[Bax_0, c0_scaling(), c0_to_c1_k(), c1_scaling(),
                        c1_to_c2_k(), c2_scaling()])
    return s.yexpr['NBD']
fitting.fit(model_func,
            [c0_scaling, c0_to_c1_k, c1_scaling, c1_to_c2_k, c2_scaling],
            ny, nt)

plt.figure()
plt.plot(nt, ny)
plt.plot(nt, model_func(nt))

print "Lag phase from dye release: %f" % params[2]
print "Main phase from dye release: %f" % params[0]
print "c0_to_c1_k, NBD: %f" % c0_to_c1_k()
print "c1_to_c2_k, NBD: %f" % c1_to_c2_k()
