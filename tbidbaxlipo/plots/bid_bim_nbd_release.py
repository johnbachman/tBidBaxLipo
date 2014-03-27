from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import fitting
from tbidbaxlipo.models.nbd.multiconf import Builder
from pysb.integrate import Solver

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

plt.figure()
params_dict = {'c1_to_c2_k': 1e-4, 'c1_scaling': 2,
               'c0_to_c1_k': 2e-3}
plt.ion()

for i in range(1, 4):
    rt = df[('Bid', 'Release', '68', i, 'TIME')].values
    ry = df[('Bid', 'Release', '68', i, 'VALUE')].values
    nt = df[('Bid', 'NBD', '68', i, 'TIME')].values
    ny = df[('Bid', 'NBD', '68', i, 'VALUE')].values

    """
    import titration_fits as tf


    plt.figure()
    plt.plot(rt, ry)

    twoexp = tf.TwoExp()
    params = twoexp.fit_timecourse(rt, ry)

    plt.plot(rt, twoexp.fit_func(rt, params))
    """
    builder = Builder(params_dict=params_dict)
    builder.build_model_multiconf(3, ny[0], normalized_data=True)

    k1_index = builder.model.parameters.index(
                        builder.model.parameters['c0_to_c1_k'])
    k2_index = builder.model.parameters.index(
                        builder.model.parameters['c1_to_c2_k'])

    #builder.model.parameters['c1_to_c2_k'].value = 1e-4
    #builder.model.parameters['c1_scaling'].value = 2.
    #builder.model.parameters['c0_to_c1_k'].value = 2e-3

    pysb_fit = fitting.fit_pysb_builder(builder, 'NBD', nt, ny)
    #print params[k1_index]
    #print params[k2_index]
    plt.plot(nt, ny)
    plt.plot(nt, pysb_fit.ypred)
    cov_x = pysb_fit.result[1]
    sd = np.sqrt(cov_x[k1_index, k1_index] * np.var(pysb_fit.residuals))
    print sd

plt.figure()
plt.hist(pysb_fit.residuals)

"""
print "Lag phase from dye release: %f" % params[2]
print "Main phase from dye release: %f" % params[0]
print "c0_to_c1_k, NBD: %f" % c0_to_c1_k()
print "c1_to_c2_k, NBD: %f" % c1_to_c2_k()
"""
