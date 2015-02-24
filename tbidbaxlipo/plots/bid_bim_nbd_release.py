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

class FitResult(object):
    """A helper class for returning and displaying fit results.

    param_dict is a dict mapping parameter name to a tuple of the
    mean and standard error for the parameter from the fit.
    """

    def __init__(self, builder, activator, nbd_site, rep_index,
                 measurement, param_dict, t, y):
        self.builder = builder
        self.activator = activator
        self.nbd_site = nbd_site
        self.rep_index = rep_index
        self.measurement = measurement
        self.param_dict = param_dict
        self.t = t
        self.y = y

    def param_result_as_string_list(self, param_name):
        """Useful for making CSV files."""
        p_mean = self.param_dict[param_name][0]
        p_se = self.param_dict[param_name][1]
        return [self.activator, self.nbd_site, str(self.rep_index),
                self.measurement, param_name, str(p_mean),
                str(p_se)]

    def timecourse_filename(self):
        return '%s_%s_r%s_%s.csv' % (self.activator, self.nbd_site, \
                                    str(self.rep_index), self.measurement)

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

def plot_endpoints(last_n_pts=3):
    nbd_sites = ['3', '5', '15', '36', '47', '54', '62', '68', '79',
                 '120', '122', '126', '138', '151', '175', '179', '184', '188']
    replicates = range(1, 4)
    activators = ['Bid', 'Bim']
    # Matrix for storing the endpoints for all mutants, replicates
    r_endpts = np.zeros((len(nbd_sites), len(replicates)))
    n_endpts = np.zeros((len(nbd_sites), len(replicates)))
    # Figure setup
    set_fig_params_for_publication()
    plt.figure('Release endpt', figsize=(4, 1.5), dpi=300) # Release figure
    plt.ylabel(r'\% Dye Release' + '\n(normalized to WT)',
               fontsize=fontsize, multialignment='center')

    plt.figure('NBD endpt', figsize=(4, 1.5), dpi=300) # NBD figure
    plt.ylabel('NBD F/$F_0$', fontsize=fontsize)
    bar_colors = {'Bid': 'gray', 'Bim': 'black'}

    # For both activators...
    for act_ix, activator in enumerate(activators):
        # Get the wild type release as a baseline, averaging over last n pts
        wt_endpts = []
        for rep_num in replicates:
            ry = df[(activator, 'Release', 'WT', rep_num, 'VALUE')].values
            wt_endpts.append(np.mean(ry[-last_n_pts:]))
        wt_mean = np.mean(wt_endpts)
        wt_sd = np.std(wt_endpts, ddof=1)

        # Now iterate over all of the mutants
        for nbd_index, nbd_site in enumerate(nbd_sites):
            # Iterate over the replicates for this mutant
            # Note that rep_num is the 1-indexed number of the replicate
            # (1, 2, 3) whereas the index is the 0-based index into the array
            # for the replicates (0, 1, 2)
            for rep_index, rep_num in enumerate(replicates):
                # Get the release data
                ry = df[(activator, 'Release', nbd_site,
                        rep_num, 'VALUE')].values
                ny = df[(activator, 'NBD', nbd_site, rep_num, 'VALUE')].values
                # Fill in entry for this replicate with mean over last n pts
                r_endpts[nbd_index, rep_index] = np.mean(ry[-last_n_pts:])
                n_endpts[nbd_index, rep_index] = np.mean(ny[-last_n_pts:])

            # Bar plot of release endpoint
            plt.figure('Release endpt')
            # Calculate percent release relative to wild type
            rel_mean = np.mean(r_endpts[nbd_index, :])
            rel_sd = np.std(r_endpts[nbd_index, :], ddof=1)
            (rel_norm_mean, rel_norm_sd) = \
                        calc_ratio_mean_sd(rel_mean, rel_sd, wt_mean, wt_sd)
            plt.bar(range(nbd_index*3 + act_ix, (nbd_index*3) + 1 + act_ix),
                    rel_norm_mean * 100,
                    width=1, color=bar_colors[activator], linewidth=0,
                    ecolor='k', capsize=1.5, yerr=rel_norm_sd * 100)

            # Bar plot of NBD endpoint
            plt.figure('NBD endpt')
            plt.bar(range(nbd_index*3 + act_ix, (nbd_index*3) + 1 + act_ix),
                    np.mean(n_endpts[nbd_index, :]),
                    width=1, color=bar_colors[activator], linewidth=0,
                    ecolor='k', capsize=1.5,
                    yerr=np.std(n_endpts[nbd_index, :], ddof=1))

    fig_names = ['Release endpt', 'NBD endpt']
    for fig_name in fig_names:
        plt.figure(fig_name)
        plt.subplots_adjust(left=0.11, bottom=0.10, right=0.97, top=0.94)
        ax = plt.gca()
        ax.set_xlim([0, len(nbd_sites) * 3])
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_tick_params(which='both', labelsize=fontsize, pad=2,
                length=2, width=0.5)
        ax.yaxis.set_tick_params(which='both', direction='out',
                        labelsize=fontsize, pad=0, length=2, width=0.5)
        ax.xaxis.labelpad = 2
        ax.yaxis.labelpad = 2
        ax.set_xticks(np.arange(1, 1 + len(nbd_sites) * 3, 3))
        ax.set_xticklabels(nbd_sites)

def plot_initial_rates(activator, num_pts=4):
    #nbd_sites = ['3', '5', '15']
    nbd_sites = ['WT', '3', '5', '15', '36', '47', '54', '62', '68', '79',
                 '120', '122', '126', '138', '151', '175', '179', '184', '188']
    replicates = range(1, 4)
    count = 0
    fit_results =[]

    for nbd_index, nbd_site in enumerate(nbd_sites):
        rn_ratios = []
        nr_ratios = []
        r_maxs = []
        n_maxs = []
        rn_errs = []
        nr_errs = []
        r_errs = []
        n_errs = []

        for rep_index in replicates:
            rt = df[(activator, 'Release', nbd_site, rep_index, 'TIME')].values
            ry = df[(activator, 'Release', nbd_site, rep_index, 'VALUE')].values
            r_lin = scipy.stats.linregress(rt[0:num_pts], ry[0:num_pts])
            r_int = r_lin[1]
            r_max = np.max(ry)
            r_slope = r_lin[0]
            r_maxs.append(r_slope)
            r_slope_err = r_lin[4]
            r_errs.append(r_slope_err)

            tb_param_dict = {'Initial rate (first %d pts)' % num_pts :
                             (r_slope, r_slope_err)}
            fit_results.append(FitResult(None, activator, nbd_site,
                                        rep_index, 'Tb release', tb_param_dict,
                                        None, None))

            if nbd_site == 'WT':
                n_maxs.append(0)
                n_errs.append(0)
                rn_ratios.append(0)
                rn_errs.append(0)
                nr_ratios.append(0)
                nr_errs.append(0)
            else:
                nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
                ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values
                # Fit line to first 10 pts
                n_lin = scipy.stats.linregress(nt[0:num_pts], ny[0:num_pts])
                n_max = np.max(ny)
                n_int = n_lin[1]
                n_slope = n_lin[0]
                n_maxs.append(n_slope)
                n_slope_err = n_lin[4]
                n_errs.append(n_slope_err)

                rn_ratio = r_slope / n_slope
                rn_err = calc_ratio_sd(r_slope, r_slope_err, n_slope,
                                       n_slope_err)
                rn_ratios.append(rn_ratio)
                rn_errs.append(rn_err)

                nr_ratio = n_slope / r_slope
                nr_err = calc_ratio_sd(n_slope, n_slope_err, r_slope,
                                       r_slope_err)
                nr_ratios.append(nr_ratio)
                nr_errs.append(nr_err)

                nbd_param_dict = {'Initial rate (first %d pts)' % num_pts :
                                  (n_slope, n_slope_err)}
                fit_results.append(FitResult(None, activator, nbd_site,
                                        rep_index, 'NBD', nbd_param_dict,
                                        None, None))

            #print "%s, rep %d, Tb slope: %f" % (nbd_site, rep_index, r_slope)
            #print "%s, rep %d, NBD slope: %f" % (nbd_site, rep_index, n_slope)

            plt.figure('%s, NBD-%s-Bax initial rates' % (activator, nbd_site),
                       figsize=(12, 5))
            plt.subplot(1, 2, 1)
            plt.plot(rt[0:num_pts], ry[0:num_pts], linestyle='', marker='.',
                     color=rep_colors[rep_index])
            plt.plot(rt[0:num_pts], r_int + r_slope * rt[0:num_pts],
                     color=rep_colors[rep_index],
                     label='%s Rep %d' % (activator, rep_index))
            plt.title('NBD-%s-Bax, Tb initial rate' % (nbd_site))
            plt.legend(loc='lower right')

            if nbd_site != 'WT':
                plt.subplot(1, 2, 2)
                plt.plot(nt[0:num_pts], ny[0:num_pts], linestyle='', marker='.',
                         color=rep_colors[rep_index])
                plt.plot(nt[0:num_pts], n_int + n_slope * nt[0:num_pts],
                         color=rep_colors[rep_index],
                         label='%s Rep %d' % (activator, rep_index))
                plt.xlabel('Time (sec)')
                plt.ylabel('$F/F_0$')
                plt.title('NBD-%s-Bax, NBD initial rate' % (nbd_site))
                plt.legend(loc='lower right')

        # NBD/Tb ratio
        plt.figure("NBD/Tb initial rate ratio", figsize=(12,6))
        plt.title("NBD/Tb initial rate ratio")
        plt.ylim((-0.08, 0.58))
        if activator == 'Bid':
            plt.bar(range(nbd_index*7, (nbd_index*7) + 3), nr_ratios,
                    width=1, color='r', ecolor='k',
                    yerr=nr_errs)
        else:
            plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), nr_ratios,
                    width=1, color='g', ecolor='k',
                    yerr=nr_errs)
        # Tb/NBD ratio
        plt.figure("Tb/NBD initial rate ratio", figsize=(12,6))
        plt.title("Tb/NBD initial rate ratio")
        plt.ylim((-100, 120))
        if activator == 'Bid':
            plt.bar(range(nbd_index*7, (nbd_index*7) + 3), rn_ratios,
                    width=1, color='r', ecolor='k',
                    yerr=rn_errs)
        else:
            plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), rn_ratios,
                    width=1, color='g', ecolor='k',
                    yerr=rn_errs)
        # Tb initial rate
        plt.figure("Tb initial rate", figsize=(12, 6))
        plt.title('Tb initial rate')
        if activator == 'Bid':
            plt.bar(range(nbd_index*7, (nbd_index*7) + 3), r_maxs,
                    width=1, color='r', ecolor='k',
                    yerr=r_errs)
        else:
            plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), r_maxs,
                    width=1, color='g', ecolor='k',
                    yerr=r_errs)
        # NBD initial rate
        plt.figure("NBD initial rate", figsize=(12, 6))
        plt.title('NBD initial rate')
        if activator == 'Bid':
            plt.bar(range(nbd_index*7, (nbd_index*7) + 3), n_maxs,
                    width=1, color='r', ecolor='k',
                    yerr=n_errs)
        else:
            plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), n_maxs,
                    width=1, color='g', ecolor='k',
                    yerr=n_errs)

    num_sites = len(nbd_sites)
    fig_names = ["NBD/Tb initial rate ratio", "Tb/NBD initial rate ratio",
                 "Tb initial rate", "NBD initial rate"]
    for fig_name in fig_names:
        plt.figure(fig_name)
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
    plt.ion()
    plot_endpoints()
    #plot_initial_rates('Bid')
    #plot_initial_rates('Bim')
    #plot_filtered_data('Bid')
    #fr = plot_3conf_fits('Bid')
    #plot_3conf_fits('Bid')
    #plot_peak_slopes('Bim')
    #plot_derivatives('Bim')
