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
import scipy.stats
from tbidbaxlipo.util.error_propagation import calc_ratio_sd

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
    nbd_sites = ['122']
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

            plt.figure()
            plt.plot(rt, ry)

            plt.figure()
            plt.plot(nt, ny)

            plt.figure()
            plt.plot(nt, ry / ny)

            ry_norm = (ry - np.min(ry)) / (np.max(ry) - np.min(ry))
            ny_norm = (ny - np.min(ny)) / (np.max(ny) - np.min(ny))

            plt.figure()
            plt.plot(rt, ry_norm, color='g')
            plt.plot(rt, ny_norm, color='b')

            plt.figure()
            plt.plot(rt, ry)
            twoexp = tf.TwoExpLinear()
            #twoexp = tf.TwoExp()
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

def plot_initial_rates(activator):
    plt.ion()
    #nbd_sites = ['15', '79']
    nbd_sites = ['3', '5', '15', '36', '47', '54', '62', '68', '79', '120',
                 '122', '126', '138', '151', '175', '179', '184', '188']
    replicates = range(1, 4)
    count = 0
    num_pts = 4
    for nbd_index, nbd_site in enumerate(nbd_sites):
        rn_ratios = []
        rn_errs = []
        for rep_index in replicates:
            rt = df[(activator, 'Release', nbd_site, rep_index, 'TIME')].values
            ry = df[(activator, 'Release', nbd_site, rep_index, 'VALUE')].values
            nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
            ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values
            # Fit line to first 10 pts
            r_lin = scipy.stats.linregress(rt[0:num_pts], ry[0:num_pts])
            n_lin = scipy.stats.linregress(nt[0:num_pts], ny[0:num_pts])
            r_int = r_lin[1]
            r_max = np.max(ry)
            n_max = np.max(ny)
            r_slope = r_lin[0] / r_max
            n_int = n_lin[1]
            n_slope = n_lin[0] / n_max

            print "%s, rep %Vd, Tb slope: %f" % (nbd_site, rep_index, r_slope)
            print "%s, rep %d, NBD slope: %f" % (nbd_site, rep_index, n_slope)
            rn_ratio = r_slope / n_slope
            rn_err = calc_ratio_sd(r_slope, r_lin[4], n_slope, n_lin[4])
            rn_ratios.append(rn_ratio)
            rn_errs.append(rn_err)

            """
            plt.figure()
            plt.plot(nt[0:num_pts], ny[0:num_pts])
            plt.plot(nt[0:num_pts], n_int + n_slope * nt[0:num_pts])
            plt.title('%s, rep %d, NBD initial rate' % (nbd_site, rep_index))

            plt.figure()
            plt.plot(rt[0:num_pts], ry[0:num_pts])
            plt.plot(rt[0:num_pts], r_int + r_slope * rt[0:num_pts])
            plt.title('%s, rep %d, Tb initial rate' % (nbd_site, rep_index))
            """

        plt.figure("Tb/NBD ratio")
        if activator == 'Bid':
            plt.bar(range(nbd_index*7, (nbd_index*7) + 3), rn_ratios,
                    width=1, color='r')
                    #yerr=rn_errs,
        else:
            plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), rn_ratios,
                    width=1, color='g')
                    #yerr=rn_errs,

    num_sites = len(nbd_sites)
    plt.figure("Tb/NBD ratio")
    ax = plt.gca()
    ax.set_xticks(np.arange(3, 3 + num_sites * 7, 7))
    ax.set_xticklabels(nbd_sites)

def plot_derivatives(activator):
    plt.ion()
    #nbd_sites = ['15', '79']
    nbd_sites = ['WT', '3', '5', '15', '36', '40', '47', '54', '62', '68',
                 '79', '120', '122', '126', '138', '151', '175', '179',
                 '184', '188']
    replicates = range(1, 4)
    count = 0
    num_pts = 4
    window = 5 # For moving average
    for nbd_index, nbd_site in enumerate(nbd_sites):
        rn_ratios = []
        r_maxs = []
        n_maxs = []
        #rn_errs = []
        for rep_index in replicates:

            rt = df[(activator, 'Release', nbd_site, rep_index, 'TIME')].values
            ry = df[(activator, 'Release', nbd_site, rep_index, 'VALUE')].values
            r_diff = np.diff(ry)
            r_avg = moving_average(r_diff, n=window)
            r_max = np.max(np.abs(r_avg))
            r_maxs.append(r_max)

            if nbd_site == 'WT':
                n_maxs.append(0)
                rn_ratios.append(0)
            else:
                nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
                ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values
                n_diff = np.diff(ny)
                n_avg = moving_average(n_diff, n=window)
                n_max = np.max(np.abs(n_avg))
                n_maxs.append(n_max)

                rn_ratio = np.max(r_avg) / np.max(n_avg)
                rn_ratios.append(rn_ratio)
            #rn_err = calc_ratio_sd(r_slope, r_lin[4], n_slope, n_lin[4])
            #rn_errs.append(rn_err)

            """
            plt.figure()
            plt.plot(nt[1:], n_diff, linestyle='', marker='.')
            plt.plot(nt[1+window-1:], n_avg)
            plt.xlabel('Time (sec)')
            plt.ylabel('dNBD/dt (F/F0 sec^-1)')
            plt.title('%s, rep %d, NBD derivative' % (nbd_site, rep_index))

            plt.figure()
            plt.plot(rt[1:], r_diff, linestyle='', marker='.')
            plt.plot(rt[1+window-1:], r_avg)
            plt.ylabel('dRel/dt (%rel sec^-1)')
            plt.title('%s, rep %d, Tb derivative' % (nbd_site, rep_index))
            """
        plt.figure("Tb/NBD peak slope ratio")
        if activator == 'Bid':
            plt.bar(range(nbd_index*7, (nbd_index*7) + 3), rn_ratios,
                    width=1, color='r')
                    #yerr=rn_errs,
        else:
            plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), rn_ratios,
                    width=1, color='g')
                    #yerr=rn_errs,
        plt.figure("Tb peak slope")
        if activator == 'Bid':
            plt.bar(range(nbd_index*7, (nbd_index*7) + 3), r_maxs,
                    width=1, color='r')
        else:
            plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), r_maxs,
                    width=1, color='g')
        plt.figure("NBD peak slope")
        if activator == 'Bid':
            plt.bar(range(nbd_index*7, (nbd_index*7) + 3), n_maxs,
                    width=1, color='r')
        else:
            plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), n_maxs,
                    width=1, color='g')

    num_sites = len(nbd_sites)
    fig_names = ["Tb/NBD peak slope ratio", "Tb peak slope", "NBD peak slope"]
    for fig_name in fig_names:
        plt.figure(fig_name)
        ax = plt.gca()
        ax.set_xticks(np.arange(3, 3 + num_sites * 7, 7))
        ax.set_xticklabels(nbd_sites)


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

if __name__ == '__main__':
    #plot_3conf_fits()
    #plot_initial_rates('Bid')
    #plot_initial_rates('Bim')
    plot_derivatives('Bid')
    plot_derivatives('Bim')
