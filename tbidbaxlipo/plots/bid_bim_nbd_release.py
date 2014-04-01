from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import fitting
from tbidbaxlipo.models.nbd.multiconf import Builder
from pysb.integrate import Solver
import itertools
import titration_fits as tf

font = {'size': 8}
matplotlib.rc('font', **font)

line_colors = {'Bid': 'r', 'Bim': 'b'}
line_styles = {1:':', 2:'-', 3:'--'}
rep_colos = {1:'r', 2:'g', 3:'b'}

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

params_dict = {'c1_to_c2_k': 1e-4, 'c1_scaling': 2,
               'c0_to_c1_k': 2e-3}

def plot_3conf_fits():
    plt.ion()
    nbd_sites = ['15']
    #nbd_sites = ['3', '5', '15', '36', '47', '54', '62', '68', '79', '120',
    #             '122', '126', '138', '151', '175', '179', '184', '188']
    replicates = range(1, 4)
    count = 0

    plt.figure("Fitted K1")
    for nbd_index, nbd_site in enumerate(nbd_sites):
        k1_means = []
        k2_means = []
        k1_sds = []
        k2_sds = []
        for rep_index in replicates:
            rt = df[('Bid', 'Release', nbd_site, rep_index, 'TIME')].values
            ry = df[('Bid', 'Release', nbd_site, rep_index, 'VALUE')].values
            nt = df[('Bid', 'NBD', nbd_site, rep_index, 'TIME')].values
            ny = df[('Bid', 'NBD', nbd_site, rep_index, 'VALUE')].values

            #deriv = np.diff(ry)
            #plt.figure()
            #plt.plot(rt[1:], deriv)
            #print rt[np.argmax(deriv)]
            #continue

            plt.figure()
            plt.plot(rt, ry)
            #twoexp = tf.TwoExpLinear()
            twoexp = tf.TwoExp()
            params = twoexp.fit_timecourse(rt, ry)
            plt.plot(rt, twoexp.fit_func(rt, params))
            print params[2]

            builder = Builder(params_dict=params_dict)
            builder.build_model_multiconf(3, ny[0], normalized_data=True)

            k1 = builder.model.parameters['c0_to_c1_k']
            k2 = builder.model.parameters['c1_to_c2_k']
            k1_index = builder.model.parameters.index(k1)
            k2_index = builder.model.parameters.index(k2)
            k1_est_index = builder.estimate_params.index(k1)
            k2_est_index = builder.estimate_params.index(k2)

            pysb_fit = fitting.fit_pysb_builder(builder, 'NBD', nt, ny)
            plt.figure()
            plt.plot(nt, ny)
            plt.plot(nt, pysb_fit.ypred)
            cov_x = pysb_fit.result[1]
            # Calculate stderr of parameters (doesn't account for covariance)
            k1_means.append(pysb_fit.params[k1_index])
            k2_means.append(pysb_fit.params[k2_index])
            k1_sds.append(np.sqrt(cov_x[k1_est_index, k1_est_index] *
                          np.var(pysb_fit.residuals)))
            k2_sds.append(np.sqrt(cov_x[k2_est_index, k2_est_index] *
                          np.var(pysb_fit.residuals)))

            plt.figure()
            s = Solver(builder.model, nt)
            s.run(param_values=pysb_fit.params)
            plt.plot(nt, s.yobs['Bax_c0'])
            plt.plot(nt, s.yobs['Bax_c1'])
            plt.plot(nt, s.yobs['Bax_c2'])

            count += 1
        plt.title(nbd_site)
        plt.xlabel('Time (sec)')
        plt.ylabel('$F/F_0$')

        plt.figure("Fitted K1")
        plt.bar(range(nbd_index*7, (nbd_index*7) + 3), k1_means,
                yerr=k1_sds, width=1, color='r')
        plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), k2_means,
                yerr=k2_sds, width=1, color='g')

    num_sites = len(nbd_sites)
    plt.figure("Fitted K1")
    ax = plt.gca()
    ax.set_xticks(np.arange(3, 3 + num_sites * 7, 7))
    ax.set_xticklabels(nbd_sites)
