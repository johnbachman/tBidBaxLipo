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
import scipy.signal
from tbidbaxlipo.util.error_propagation import calc_ratio_sd, calc_ratio_mean_sd
from tbidbaxlipo.util import moving_average
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis, \
                             fontsize

font = {'size': 10}
matplotlib.rc('font', **font)

line_colors = {'Bid': 'r', 'Bim': 'b'}
line_styles = {1:':', 2:'-', 3:'--'}
rep_colors = {1:'r', 2:'g', 3:'b'}

def _mean_sd(p_name, builder, pysb_fit):
    """Given the named parameter and the fit result, get mean and SE."""
    p = builder.model.parameters[p_name]
    p_index = builder.model.parameters.index(p)
    p_est_index = builder.estimate_params.index(p)
    p_mean = pysb_fit.params[p_index]
    cov_x = pysb_fit.result[1]
    if cov_x is not None:
        p_sd = np.sqrt(cov_x[p_est_index, p_est_index] *
                        np.var(pysb_fit.residuals))
    else:
        p_sd = np.nan

    return (p_mean, p_sd)


params_dict = {'c1_to_c2_k': 1e-4, 'c1_scaling': 2,
               'c0_to_c1_k': 2e-3}

def plot_2conf_fits(activator):
    #nbd_sites = ['15', '54', '62', '68', '79', '126', '138', '175']
    nbd_sites = ['3', '5', '15', '36', '47', '54', '62', '68', '79', '120',
                 '122', '126', '138', '151', '175', '179', '184', '188']
    replicates = range(1, 4)
    count = 0
    fit_results = []

    plt.figure("Fitted k1")
    for nbd_index, nbd_site in enumerate(nbd_sites):
        k1_means = []
        k1_sds = []

        for rep_index in replicates:
            nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
            ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values

            plt.figure('%s, NBD-%s-Bax Fits' % (activator, nbd_site),
                       figsize=(6, 5))

            builder = Builder(params_dict=params_dict)
            builder.build_model_multiconf(2, ny[0], normalized_data=True)
            # Rough guesses for parameters
            # Guess that the final scaling is close to the final data value
            builder.model.parameters['c1_scaling'].value = ny[-1]
            # Rough guess for the timescale
            builder.model.parameters['c0_to_c1_k'].value = 1e-3

            pysb_fit = fitting.fit_pysb_builder(builder, 'NBD', nt, ny)
            plt.plot(nt, ny, linestyle='', marker='.',
                    color=rep_colors[rep_index])
            plt.plot(nt, pysb_fit.ypred,
                     label='%s Rep %d' % (activator, rep_index),
                     color=rep_colors[rep_index])
            plt.xlabel('Time (sec)')
            plt.ylabel('$F/F_0$')
            plt.title('NBD $F/F_0$ fits, NBD-%s-Bax' % nbd_site)
            plt.legend(loc='lower right')
            # Calculate stderr of parameters (doesn't account for covariance)
            (k1_mean, k1_sd) = _mean_sd('c0_to_c1_k', builder, pysb_fit)
            k1_means.append(k1_mean)
            k1_sds.append(k1_sd)
            (c1_mean, c1_sd) = _mean_sd('c1_scaling', builder, pysb_fit)

            param_dict = {'c0_to_c1_k': (k1_mean, k1_sd),
                          'c1_scaling': (c1_mean, c1_sd)}
            fit_results.append(FitResult(builder, activator, nbd_site,
                                         rep_index, 'NBD', param_dict,
                                         nt, pysb_fit.ypred))

            count += 1

        plt.figure('%s, NBD-%s-Bax Fits' % (activator, nbd_site))
        plt.tight_layout()

        plt.figure("Fitted k1")
        plt.bar(range(nbd_index*4, (nbd_index*4) + 3), k1_means,
                yerr=k1_sds, width=1, color='r', ecolor='k')

    num_sites = len(nbd_sites)
    plt.figure("Fitted k1")
    plt.ylabel('Fitted k1 ($sec^{-1}$)')
    ax = plt.gca()
    ax.set_xticks(np.arange(1.5, 1.5 + num_sites * 4, 4))
    ax.set_xticklabels(nbd_sites)

    return fit_results

def plot_3conf_fits(activator):
    #nbd_sites = ['15', '54', '62', '68', '79', '126', '138', '175']
    nbd_sites = ['3', '5', '15', '36', '47', '54', '62', '68', '79', '120',
                 '122', '126', '138', '151', '175', '179', '184', '188']
    replicates = range(1, 4)
    count = 0
    fit_results = []
    plt.figure("Fitted k1/k2")
    for nbd_index, nbd_site in enumerate(nbd_sites):
        k1_means = []
        k2_means = []
        k1_sds = []
        k2_sds = []

        for rep_index in replicates:
            rt = df[(activator, 'Release', nbd_site, rep_index, 'TIME')].values
            ry = df[(activator, 'Release', nbd_site, rep_index, 'VALUE')].values
            nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
            ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values

            """
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
            """
            plt.figure('%s, NBD-%s-Bax Fits' % (activator, nbd_site),
                       figsize=(12, 5))
            plt.subplot(1, 2, 1)
            plt.plot(rt, ry,
                     linestyle='', marker='.',
                     color=rep_colors[rep_index])
            twoexp = tf.TwoExpLinear()
            #twoexp = tf.TwoExp()
            params = twoexp.fit_timecourse(rt, ry)
            plt.plot(rt, twoexp.fit_func(rt, params),
                     label='%s Rep %d' % (activator, rep_index),
                     color=rep_colors[rep_index])
            plt.xlabel('Time (sec)')
            plt.ylabel('% Tb release')
            plt.title('%% Tb release fits, NBD-%s-Bax' % nbd_site)
            plt.legend(loc='lower right')

            builder = Builder(params_dict=params_dict)
            builder.build_model_multiconf(3, ny[0], normalized_data=True)
            # Add some initial guesses
            # Guess that the scaling value for the final conformation is close
            # to the final data value
            builder.model.parameters['c2_scaling'].value = ny[-1]
            # Guess that the scaling value for the intermediate conformation
            # is close to the value at ~300 sec
            c1_timescale_seconds = 300
            c1_timescale_index = np.where(nt > 300)[0].min()
            builder.model.parameters['c1_scaling'].value = \
                                                    ny[c1_timescale_index]
            # Rough guesses for the timescales of the first and second
            # transitions
            builder.model.parameters['c0_to_c1_k'].value = 5e-3
            builder.model.parameters['c1_to_c2_k'].value = 5e-4

            pysb_fit = fitting.fit_pysb_builder(builder, 'NBD', nt, ny)
            plt.subplot(1, 2, 2)
            plt.plot(nt, ny, linestyle='', marker='.',
                    color=rep_colors[rep_index])
            plt.plot(nt, pysb_fit.ypred,
                     label='%s Rep %d' % (activator, rep_index),
                     color=rep_colors[rep_index])
            plt.xlabel('Time (sec)')
            plt.ylabel('$F/F_0$')
            plt.title('NBD $F/F_0$ fits, NBD-%s-Bax' % nbd_site)
            plt.legend(loc='lower right')
            # Calculate stderr of parameters (doesn't account for covariance)
            (k1_mean, k1_sd) = _mean_sd('c0_to_c1_k', builder, pysb_fit)
            (k2_mean, k2_sd) = _mean_sd('c1_to_c2_k', builder, pysb_fit)
            (c1_mean, c1_sd) = _mean_sd('c1_scaling', builder, pysb_fit)
            (c2_mean, c2_sd) = _mean_sd('c2_scaling', builder, pysb_fit)

            param_dict = {'c0_to_c1_k': (k1_mean, k1_sd),
                          'c1_scaling': (c1_mean, c1_sd),
                          'c1_to_c2_k': (k2_mean, k2_sd),
                          'c2_scaling': (c2_mean, c2_sd)}

            fit_results.append(FitResult(builder, activator, nbd_site,
                                         rep_index, 'NBD', param_dict,
                                         nt, pysb_fit.ypred))

            k1_means.append(k1_mean)
            k2_means.append(k2_mean)
            k1_sds.append(k1_sd)
            k2_sds.append(k2_sd)

            """
            plt.figure()
            s = Solver(builder.model, nt)
            s.run(param_values=pysb_fit.params)
            plt.plot(nt, s.yobs['Bax_c0'])
            plt.plot(nt, s.yobs['Bax_c1'])
            plt.plot(nt, s.yobs['Bax_c2'])
            """
            count += 1
        #plt.title(nbd_site)
        #plt.xlabel('Time (sec)')
        #plt.ylabel('$F/F_0$')
        plt.figure('%s, NBD-%s-Bax Fits' % (activator, nbd_site))
        plt.tight_layout()

        plt.figure("Fitted k1/k2")
        plt.bar(range(nbd_index*7, (nbd_index*7) + 3), k1_means,
                yerr=k1_sds, width=1, color='r', ecolor='k')
        plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), k2_means,
                yerr=k2_sds, width=1, color='g', ecolor='k')

    num_sites = len(nbd_sites)
    plt.figure("Fitted k1/k2")
    plt.ylabel('Fitted k1/k2 ($sec^{-1}$)')
    ax = plt.gca()
    ax.set_xticks(np.arange(3, 3 + num_sites * 7, 7))
    ax.set_xticklabels(nbd_sites)

    return fit_results

def plot_peak_slopes(activator):
    #nbd_sites = ['WT', '3', '5', '54', '62', '68']
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
            r_avg = moving_average(ry, n=window)
            r_diff = np.diff(r_avg)
            r_max = np.max(np.abs(r_diff))
            r_maxs.append(r_max)

            if nbd_site == 'WT':
                n_maxs.append(0)
                rn_ratios.append(0)
            else:
                nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
                ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values
                n_avg = moving_average(ny, n=window)
                n_diff = np.diff(n_avg)
                n_max = np.max(np.abs(n_diff))
                n_maxs.append(n_max)

                rn_ratio = r_max / n_max
                rn_ratios.append(rn_ratio)
            #rn_err = calc_ratio_sd(r_slope, r_lin[4], n_slope, n_lin[4])
            #rn_errs.append(rn_err)

            plt.figure('%s, NBD-%s-Bax derivative' % (activator, nbd_site),
                       figsize=(12, 5))
            plt.subplot(1, 2, 1)
            plt.plot(rt[1+window-1:], r_diff, color=rep_colors[rep_index],
                     label='%s Rep %d' % (activator, rep_index))
            #plt.plot(rt[1+window-1:], r_avg)
            plt.ylabel('dRel/dt (% rel $sec^{-1}$)')
            plt.title('NBD-%s-Bax, Tb derivative' % nbd_site)
            plt.legend(loc='upper right')

            if nbd_site != 'WT':
                plt.subplot(1, 2, 2)
                plt.plot(nt[1+window-1:], n_diff, color=rep_colors[rep_index],
                         label='%s Rep %d' % (activator, rep_index))
                #plt.plot(nt[1+window-1:], n_avg)
                plt.xlabel('Time (sec)')
                plt.ylabel('dNBD/dt ($F/F_0\ sec^{-1}$)')
                plt.title('NBD-%s-Bax, NBD derivative' % nbd_site)
                plt.legend(loc='upper right')

                # Plot normalized slopes
                plt.figure('%s, NBD-%s-Bax normalized derivative' %
                           (activator, nbd_site))
                n_diff_norm = n_diff / np.max(np.abs(n_diff))
                r_diff_norm = r_diff / np.max(np.abs(r_diff))
                plt.plot(rt[1+window-1:], r_diff_norm, color=rep_colors[rep_index],
                        linestyle=line_styles[2],
                        label='%s Rep %d' % (activator, rep_index))
                plt.plot(nt[1+window-1:], n_diff_norm, color=rep_colors[rep_index],
                        linestyle=line_styles[3])
                plt.xlabel('Time (sec)')
                plt.ylabel('% max rate')
                plt.title('%s, NBD-%s-Bax normalized derivative' %
                          (activator, nbd_site))
                plt.legend(loc='upper right')

        plt.figure('%s, NBD-%s-Bax derivative' % (activator, nbd_site))
        plt.tight_layout()
        plt.show()
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
        """
    """
    num_sites = len(nbd_sites)
    fig_names = ["Tb/NBD peak slope ratio", "Tb peak slope", "NBD peak slope"]
    for fig_name in fig_names:
        plt.figure(fig_name)
        ax = plt.gca()
        ax.set_xticks(np.arange(3, 3 + num_sites * 7, 7))
        ax.set_xticklabels(nbd_sites)
    """

def plot_filtered_data(activator):
    nbd_sites = ['15', '54']
    #nbd_sites = ['WT', '3', '5', '15', '36', '40', '47', '54', '62', '68',
    #             '79', '120', '122', '126', '138', '151', '175', '179',
    #             '184', '188']
    replicates = range(1, 4)
    count = 0
    num_pts = 4
    window = 5 # For moving average
    # Create an order 3 lowpass butterworth filter.
    b, a = scipy.signal.butter(1, 0.2)

    for nbd_index, nbd_site in enumerate(nbd_sites):
        for rep_index in replicates:
            """
            rt = df[(activator, 'Release', nbd_site, rep_index, 'TIME')].values
            ry = df[(activator, 'Release', nbd_site, rep_index, 'VALUE')].values
            r_filt = scipy.signal.filtfilt(b, a, ry)
            plt.figure()
            plt.plot(rt, ry, linestyle='', marker='.')
            plt.plot(rt, r_filt)
            """
            nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
            ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values
            n_filt = scipy.signal.filtfilt(b, a, ny)
            plt.figure()
            plt.plot(nt, ny, linestyle='', marker='.')
            plt.plot(nt, n_filt)

if __name__ == '__main__':
    #nbd_sites = ['3', '5', '15', '36', '47', '54', '62', '68', '79',
    #             '120', '122', '126', '138', '151', '175', '179', '184', '188']
    #plot_all(df, nbd_residues, file_basename='data1_raw')
    plot_endpoints(df, nbd_residues, file_basename='data1')
    #plot_initial_rates('Bid')
    #plot_initial_rates('Bim')
    #plot_filtered_data('Bid')
    #fr = plot_3conf_fits('Bid')
    #plot_3conf_fits('Bid')
    #plot_peak_slopes('Bim')
    #plot_derivatives('Bim')
