import csv
import itertools
import collections
from matplotlib import pyplot as plt
import numpy as np
import scipy.stats
import scipy.signal # For low-pass filtering in derivatives function
from tbidbaxlipo.util import set_fig_params_for_publication, fontsize, \
                             format_axis
from tbidbaxlipo.util.error_propagation import calc_ratio_mean_sd
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from tbidbaxlipo.util import fitting, moving_average
from tbidbaxlipo.models.nbd.multiconf import Builder
import tbidbaxlipo.plots.titration_fits as tf
import tbidbaxlipo.util.calculate_error_variance as cev

line_colors = {'Bid': 'r', 'Bim': 'b'}
dtype_line_colors = {'Release':'r', 'NBD':'g', 'FRET':'b'}
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

def plot_all(df, nbd_residues, datatypes, activators=None, replicates=None,
             file_basename=None, normalize_nbd=False):
    """Release should be first in the list of datatypes."""
    if datatypes[0] != 'Release':
        raise ValueError("Release should be first in the datatypes list.")
    # Some definitions
    # Define some default values
    if activators is None:
        activators = ['Bid', 'Bim']
    if replicates is None:
        replicates = (1, 2, 3)
    num_subplots = len(datatypes)
    # Labels and titles for each datatype
    ylabels = {'Release': '% Dye Release',
               'NBD': 'F/$F_0$',
               'FRET': '% FRET'}
    # Every mutant gets its own plot
    for nbd_index, nbd_site in enumerate(nbd_residues):
        if len(datatypes) == 2:
            fig_width = 11
        elif len(datatypes) == 3:
            fig_width = 14
        fig = plt.figure(figsize=(fig_width, 5))
        fig.set_tight_layout(True)
        # Define the subplot titles
        titles = {'Release': r'Dye release for NBD-%s-Bax' % nbd_site,
                  'NBD': 'NBD F/$F_0$ for NBD-%s-Bax' % nbd_site,
                  'FRET': 'FRET, NBD-%s-Bax' % nbd_site}
        # Every datatype gets its own subplot
        for dtype_ix, dtype in enumerate(datatypes):
            # There is no NBD/FRET curve for WT Bax, so skip
            if (dtype == 'NBD' and nbd_site == 'WT') or \
               (dtype == 'FRET' and nbd_site == 'WT'):
                continue
            plt.subplot(1, num_subplots, dtype_ix + 1)
            # Activators and replicates are separate lines on the same plot
            for activator in activators:
                # Iterate over replicates...
                for i in replicates:
                    # Get the data
                    t = df[(activator, dtype, nbd_site, i, 'TIME')]
                    v = df[(activator, dtype, nbd_site, i, 'VALUE')]
                    # TODO Do the normalization here if desired TODO
                    if normalize_nbd:
                        v = v / v[0]
                    plt.plot(t, v, label='%s Rep %d' % (activator, i),
                            color=line_colors[activator],
                            linestyle=line_styles[i])
                    plt.xlabel('Time (sec)')
                    plt.ylabel(ylabels[dtype])
                    plt.title(titles[dtype])
                    plt.legend(loc='lower right')
            # Datatype-specific formatting
            if dtype == 'Release':
                plt.ylim([0, 100])
        if file_basename:
            plt.savefig('%s_%s.pdf' % (file_basename, nbd_site))
            plt.savefig('%s_%s.png' % (file_basename, nbd_site))

def plot_all_by_replicate(df, nbd_residues, datatypes, activators=None,
                          replicates=None, file_basename=None,
                          normalize_nbd=False):
    """Release should be first in the datatypes list."""
    if datatypes[0] != 'Release':
        raise ValueError("Release should be first in the datatypes list.")
    # Define some default values
    if activators is None:
        activators = ['Bid', 'Bim']
    if replicates is None:
        replicates = (1, 2, 3)
    # Every mutant gets its own plot
    for nbd_index, nbd_site in enumerate(nbd_residues):
        for activator in activators:
            plt.figure(figsize=(14, 5))
            # Each replicate gets its own subplot
            ax2_list = []
            nbd_max = -np.inf
            nbd_min = np.inf
            for rep in replicates:
                # Make the release plot
                plt.subplot(1, 3, rep)
                ax1 = plt.gca()
                ax2 = ax1.twinx()
                ax2_list.append(ax2)
                for dtype in datatypes:
                    # Skip WT for NBD and FRET
                    if (dtype == 'NBD' and nbd_site == 'WT') or \
                       (dtype == 'FRET' and nbd_site == 'WT'):
                        continue
                    # Get the data
                    t = df[(activator, dtype, nbd_site, rep, 'TIME')]
                    v = df[(activator, dtype, nbd_site, rep, 'VALUE')]
                    # Set the axis
                    if dtype == 'NBD':
                        if normalize_nbd:
                            v = v / v[0]
                        ax = ax2
                        if np.max(v) > nbd_max:
                            nbd_max = np.max(v)
                        if np.min(v) < nbd_min:
                            nbd_min = np.min(v)
                    elif dtype == 'Release' or dtype == 'FRET':
                        ax = ax1
                    else:
                        raise ValueError("Unknown datatype: %s" % dtype)
                    # Plot the data
                    ax.plot(t, v, label='%s, %s' % (activator, dtype),
                            color=dtype_line_colors[dtype])

                # Adjust and label the figure
                ax1.set_ylim([0, 100])
                if 'Release' in datatypes and 'FRET' in datatypes:
                    ax1.set_ylabel('% FRET, % Release')
                else:
                    ax1.set_ylabel('% Release')

                ax2.set_xlabel('Time (sec)')
                ax2.set_ylabel('NBD $F/F_0$')
                ax2.set_title('NBD-%s-Bax, %s, Rep %d' %
                              (nbd_site, activator, rep))
                ax1_lines, ax1_labels = ax1.get_legend_handles_labels()
                ax2_lines, ax2_labels = ax2.get_legend_handles_labels()
                ax2.legend(ax1_lines + ax2_lines, ax1_labels + ax2_labels,
                           loc='lower right', prop={'size':8})
                plt.xticks([0, 1000, 2000, 3000, 4000])
            # Give all plots same range for NBD
            for ax2 in ax2_list:
                ax2.set_ylim([nbd_min * 0.95, nbd_max * 1.05])
            plt.subplots_adjust(wspace=0.4, left=0.06, right=0.95)
            # Output file, if desired
            if file_basename:
                plt.savefig('%s_%s_%s.pdf' %
                            (file_basename, nbd_site, activator))
                plt.savefig('%s_%s_%s.png' %
                            (file_basename, nbd_site, activator))

def calc_barplot_width(num_sites, rmarg=0.13, lmarg=0.11):
    # I found that these numbers worked for a 4 inch wide figure with the
    # given relative left and right margins and 19 sites. This allows the
    # same proportions for endpoint plots with different numbers of sites.
    abs_left = 4 * lmarg
    abs_right = 4 * rmarg
    abs_middle = 4 * (1 - lmarg - rmarg)
    abs_per_site = abs_middle / 19.
    fig_width = abs_per_site * num_sites + abs_left + abs_right
    rel_left = abs_left / float(fig_width)
    rel_right = 1 - (abs_right / float(fig_width))
    return (fig_width, rel_left, rel_right)

def plot_nbd_endpoints(df, nbd_sites, last_n_pts=3, file_basename=None,
                       normalize_nbd=False):
    replicates = range(1, 4)
    activators = ['Bid', 'Bim']
    # Filter out the WT residue from the list, if its there
    nbd_sites_no_wt = [s for s in nbd_sites if s != 'WT']
    # Matrix for storing the endpoints for all mutants, replicates
    n_endpts = np.zeros((len(nbd_sites_no_wt), len(replicates)))
    # Figure setup
    set_fig_params_for_publication()
    if normalize_nbd:
        yaxis_label = 'NBD F/$F_0$'
        (fig_width, rel_left, rel_right) = \
                            calc_barplot_width(len(nbd_sites_no_wt))
    else:
        yaxis_label = 'NBD fluorescence (RFU)'
        # If we're not normalizing the NBD values, then the additional 0s will
        # push the axis label off the left-hand side. Therefore we adjust the
        # left margin accordingly.
        (fig_width, rel_left, rel_right) = \
                            calc_barplot_width(len(nbd_sites_no_wt), lmarg=0.13)
    plt.figure(file_basename, figsize=(fig_width, 1.5), dpi=300)
    plt.ylabel(yaxis_label, fontsize=fontsize)
    bar_colors = {'Bid': 'gray', 'Bim': 'black'}

    # For both activators...
    for act_ix, activator in enumerate(activators):
        # Now iterate over all of the mutants
        for nbd_index, nbd_site in enumerate(nbd_sites_no_wt):
            # Iterate over the replicates for this mutant
            # Note that rep_num is the 1-indexed number of the replicate
            # (1, 2, 3) whereas the index is the 0-based index into the array
            # for the replicates (0, 1, 2)
            for rep_index, rep_num in enumerate(replicates):
                # Get the release data
                ny = df[(activator, 'NBD', nbd_site, rep_num, 'VALUE')].values
                # Fill in entry for this replicate with mean over last n pts
                if normalize_nbd:
                    endpt_vals = np.mean(ny[-last_n_pts:] / ny[0])
                else:
                    endpt_vals = np.mean(ny[-last_n_pts:])
                n_endpts[nbd_index, rep_index] = endpt_vals

            # Bar plot of NBD endpoint
            plt.bar(range(nbd_index*3 + act_ix, (nbd_index*3) + 1 + act_ix),
                    np.mean(n_endpts[nbd_index, :]),
                    width=1, color=bar_colors[activator], linewidth=0,
                    ecolor='k', capsize=1.5,
                    yerr=np.std(n_endpts[nbd_index, :], ddof=1))

    # Add horizontal gridlines
    x_lbound = -1
    x_ubound = len(nbd_sites_no_wt) * 3
    plt.hlines(range(0, 6), x_lbound, x_ubound, color='0.85', zorder=1,
               linewidth=0.5)
    # Format the plot
    plt.subplots_adjust(left=rel_left, bottom=0.10, right=rel_right, top=0.94)
    ax = plt.gca()
    ax.set_xlim([x_lbound, x_ubound])
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_tick_params(which='both', labelsize=fontsize, pad=2,
            length=2, width=0.5)
    ax.yaxis.set_tick_params(which='both', direction='out',
                    labelsize=fontsize, pad=0, length=2, width=0.5)
    ax.xaxis.labelpad = 2
    ax.yaxis.labelpad = 2
    ax.set_xticks(np.arange(1, 1 + len(nbd_sites_no_wt) * 3, 3))
    ax.set_xticklabels(nbd_sites_no_wt)
    # Create legend with rectangles to match colors of bars
    bid_patch = mpatches.Patch(color=bar_colors['Bid'], label='cBid')
    bim_patch = mpatches.Patch(color=bar_colors['Bim'], label='Bim')
    leg = plt.legend(handles=[bid_patch, bim_patch],
            bbox_to_anchor=(1.02, 0.5), loc='center left', borderaxespad=0.,
            prop={'size': fontsize}, handlelength=1)
    leg.draw_frame(False)
    # Output the file, if desired
    if file_basename:
        plt.savefig('%s.pdf' % file_basename)
        plt.savefig('%s.png' % file_basename)

def plot_release_endpoints(df, nbd_sites, normalized_to_wt=False,
                           last_n_pts=3, file_basename=None):
    # Check for nonsense
    if normalized_to_wt and 'WT' not in nbd_sites:
        raise ValueError("Cannot normalized to WT release if WT is not in the "
                         "dataset!")
    replicates = range(1, 4)
    activators = ['Bid', 'Bim']
    # If we're normalizing to WT, filter the WT entry from the site list
    if normalized_to_wt:
        nbd_sites_filt = [s for s in nbd_sites if s != 'WT']
    else:
        nbd_sites_filt = nbd_sites
    # Matrix for storing the endpoints for all mutants, replicates
    r_endpts = np.zeros((len(nbd_sites_filt), len(replicates)))
    # Figure setup
    set_fig_params_for_publication()
    (fig_width, rel_left, rel_right) = \
                        calc_barplot_width(len(nbd_sites_filt))
    plt.figure(file_basename, figsize=(fig_width, 1.5), dpi=300)
    # Set the yaxis label according to whether we're normalizing
    if normalized_to_wt:
        ylabel_str = r'\% Dye Release' + '\n(normalized to WT)'
    else:
        ylabel_str = r'\% Dye Release'
    plt.ylabel(ylabel_str, fontsize=fontsize, multialignment='center')
    # Bar colors for the different activators
    bar_colors = {'Bid': 'gray', 'Bim': 'black'}

    # For both activators...
    for act_ix, activator in enumerate(activators):
        # Get the wild type release as a baseline, averaging over last n pts
        if normalized_to_wt:
            wt_endpts = []
            for rep_num in replicates:
                ry = df[(activator, 'Release', 'WT', rep_num, 'VALUE')].values
                wt_endpts.append(np.mean(ry[-last_n_pts:]))
            wt_mean = np.mean(wt_endpts)
            wt_sd = np.std(wt_endpts, ddof=1)

        # Now iterate over all of the mutants
        for nbd_index, nbd_site in enumerate(nbd_sites_filt):
            # Iterate over the replicates for this mutant
            # Note that rep_num is the 1-indexed number of the replicate
            # (1, 2, 3) whereas the index is the 0-based index into the array
            # for the replicates (0, 1, 2)
            for rep_index, rep_num in enumerate(replicates):
                # Get the release data
                ry = df[(activator, 'Release', nbd_site,
                        rep_num, 'VALUE')].values
                # Fill in entry for this replicate with mean over last n pts
                r_endpts[nbd_index, rep_index] = np.mean(ry[-last_n_pts:])

            # Bar plot of release endpoint
            # Calculate percent release relative to wild type
            rel_mean = np.mean(r_endpts[nbd_index, :])
            rel_sd = np.std(r_endpts[nbd_index, :], ddof=1)
            if normalized_to_wt:
                (rel_mean, rel_sd) = \
                            calc_ratio_mean_sd(rel_mean, rel_sd, wt_mean, wt_sd)
                rel_mean *= 100
                rel_sd *= 100

            plt.bar(range(nbd_index*3 + act_ix, (nbd_index*3) + 1 + act_ix),
                    rel_mean, width=1, color=bar_colors[activator],
                    linewidth=0, ecolor='k', capsize=1.5, yerr=rel_sd)

    # Add horizontal gridlines
    x_lbound = -1
    x_ubound = len(nbd_sites_filt) * 3
    # Format the plot
    plt.subplots_adjust(left=rel_left, bottom=0.10, right=rel_right, top=0.94)
    ax = plt.gca()
    ax.set_xlim([x_lbound, x_ubound])
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_tick_params(which='both', labelsize=fontsize, pad=2,
            length=2, width=0.5)
    ax.yaxis.set_tick_params(which='both', direction='out',
                    labelsize=fontsize, pad=0, length=2, width=0.5)
    ax.xaxis.labelpad = 2
    ax.yaxis.labelpad = 2
    ax.set_xticks(np.arange(1, 1 + len(nbd_sites_filt) * 3, 3))
    ax.set_xticklabels(nbd_sites_filt)
    # Create legend with rectangles to match colors of bars
    bid_patch = mpatches.Patch(color=bar_colors['Bid'], label='cBid')
    bim_patch = mpatches.Patch(color=bar_colors['Bim'], label='Bim')
    leg = plt.legend(handles=[bid_patch, bim_patch],
            bbox_to_anchor=(1.02, 0.5), loc='center left', borderaxespad=0.,
            prop={'size': fontsize}, handlelength=1)
    leg.draw_frame(False)
    # Add hlines every 20 percent until the top of the plot is reached
    (ymin, ymax) = ax.get_ylim()
    plt.hlines(range(20, int(ymax), 20), x_lbound, x_ubound, color='0.85',
               zorder=1, linewidth=0.5)
    # Output the file, if desired
    if file_basename:
        plt.savefig('%s.pdf' % file_basename)
        plt.savefig('%s.png' % file_basename)

def plot_initial_rate_fits(df, nbd_sites, activator, num_pts=4, plot=False):
    replicates = range(1, 4)
    fit_results =[]
    # Lists for storing all of the various fits, slopes
    r_slopes_all = np.zeros((len(nbd_sites), len(replicates)))
    n_slopes_all = np.zeros((len(nbd_sites), len(replicates)))
    r_errs_all = np.zeros((len(nbd_sites), len(replicates)))
    n_errs_all = np.zeros((len(nbd_sites), len(replicates)))
    n_max_all = np.zeros((len(nbd_sites), len(replicates)))
    # Iterate over all of the mutants
    for nbd_index, nbd_site in enumerate(nbd_sites):
        # Iterate over the replicates for this mutant
        # Note that rep_num is the 1-indexed number of the replicate
        # (1, 2, 3) whereas the index is the 0-based index into the array for
        # the replicates (0, 1, 2)
        for rep_index, rep_num in enumerate(replicates):
            # Get the release data
            rt = df[(activator, 'Release', nbd_site, rep_num, 'TIME')].values
            ry = df[(activator, 'Release', nbd_site, rep_num, 'VALUE')].values
            # Fit line to first n pts
            r_lin = scipy.stats.linregress(rt[0:num_pts], ry[0:num_pts])
            # Store the slope of the line
            r_slope = r_lin[0]
            r_slopes_all[nbd_index, rep_index] = r_slope
            r_int = r_lin[1] # Intercept
            # Uncertainty associated with the slope
            r_slope_err = r_lin[4]
            r_errs_all[nbd_index, rep_index] = r_slope_err
            # Store fit result
            tb_param_dict = {'Initial rate (first %d pts)' % num_pts :
                             (r_slope, r_slope_err)}
            fit_results.append(FitResult(None, activator, nbd_site,
                                        rep_num, 'Tb release', tb_param_dict,
                                        None, None))
            # Now do the NBD slope calculation, but ignore the WT (since there
            # is no NBD label)
            if nbd_site == 'WT':
                n_max_all[nbd_index, rep_index] = 0
                n_errs_all[nbd_index, rep_index] = 0
                n_slopes_all[nbd_index, rep_index] = 0
            else:
                # Get the NBD data
                nt = df[(activator, 'NBD', nbd_site, rep_num, 'TIME')].values
                ny = df[(activator, 'NBD', nbd_site, rep_num, 'VALUE')].values
                # Fit line to first n pts
                n_lin = scipy.stats.linregress(nt[0:num_pts], ny[0:num_pts])
                # Maximum NBD F/F0
                n_max_all[nbd_index, rep_index] = np.max(ny)
                # Slope and intercept
                n_slope = n_lin[0]
                n_slopes_all[nbd_index, rep_index] = n_slope
                n_int = n_lin[1] # Intercept
                # Uncertainty associated with slope
                n_slope_err = n_lin[4]
                n_errs_all[nbd_index, rep_index] = n_slope_err
                # Store fit result
                nbd_param_dict = {'Initial rate (first %d pts)' % num_pts :
                                  (n_slope, n_slope_err)}
                fit_results.append(FitResult(None, activator, nbd_site,
                                        rep_num, 'NBD', nbd_param_dict,
                                        None, None))
            if plot:
                # Make subpanel showing linear fit for release slope
                plt.figure('%s, NBD-%s-Bax initial rates' %
                          (activator, nbd_site), figsize=(12, 5))
                # Dye release slope subplot
                plt.subplot(1, 2, 1)
                # Plot the data
                plt.plot(rt[0:num_pts], ry[0:num_pts], linestyle='',
                         marker='o', color=rep_colors[rep_num])
                # Plot the fit
                plt.plot(rt[0:num_pts], r_int + r_slope * rt[0:num_pts],
                         color=rep_colors[rep_num],
                         label='%s Rep %d' % (activator, rep_num))
                plt.title('NBD-%s-Bax, Tb initial rate' % (nbd_site))
                plt.legend(loc='lower right')

            # Make subpanel showing linear fit for NBD slope
            if plot and nbd_site != 'WT':
                plt.subplot(1, 2, 2)
                # Plot the data
                plt.plot(nt[0:num_pts], ny[0:num_pts], linestyle='',
                         marker='o', color=rep_colors[rep_num])
                # Plot the fit
                plt.plot(nt[0:num_pts], n_int + n_slope * nt[0:num_pts],
                         color=rep_colors[rep_num],
                         label='%s Rep %d' % (activator, rep_num))
                plt.xlabel('Time (sec)')
                plt.ylabel('$F/F_0$')
                plt.title('NBD-%s-Bax, NBD initial rate' % (nbd_site))
                plt.legend(loc='lower right')

    return fit_results

class InitialRateSamples(object):
    def __init__(self, r_slopes, n_slopes, n_maxes, wt_r_slopes, time):
        self.r_slopes = r_slopes
        self.n_slopes = n_slopes
        self.n_maxes = n_maxes
        self.time = time
        self.wt_r_slopes = wt_r_slopes

    def release_avg_std(self):
        r_slope_avgs = np.mean(self.r_slopes, axis=2)
        r_slope_stds = np.std(self.r_slopes, axis=2, ddof=1)
        return (r_slope_avgs, r_slope_stds)

    def release_avg_std_wt_normalized(self):
        if np.any(np.isnan(self.wt_r_slopes)):
            raise ValueError("Can't normalize to WT because there are "
                             "missing entries for initial WT release "
                             "rates (found NaNs in the matrix).")
        (wt_r_avgs, wt_r_stds) = self.wt_release_avg_std()
        (r_slope_avgs, r_slope_stds) = self.release_avg_std()
        r_norm_avgs = np.zeros(r_slope_avgs.shape)
        r_norm_stds = np.zeros(r_slope_stds.shape)

        for act_ix in range(r_slope_avgs.shape[0]):
            for nbd_index in range(r_slope_avgs.shape[1]):
                (r_norm_avgs[act_ix, nbd_index],
                 r_norm_stds[act_ix, nbd_index]) = calc_ratio_mean_sd(
                                            r_slope_avgs[act_ix, nbd_index],
                                            r_slope_stds[act_ix, nbd_index],
                                            wt_r_avgs[act_ix],
                                            wt_r_stds[act_ix])
        return (r_norm_avgs, r_norm_stds)

    def release_slopes_wt_normalized(self):
        """Get the release for each individual replicate normalized
        to the WT release at the timepoint."""
        if np.any(np.isnan(self.wt_r_slopes)):
            raise ValueError("Can't normalize to WT because there are "
                             "missing entries for initial WT release "
                             "rates (found NaNs in the matrix).")
        (wt_r_avgs, wt_r_stds) = self.wt_release_avg_std()
        r_slopes_norm = np.zeros(self.r_slopes.shape)
        for act_ix in range(r_slopes_norm.shape[0]):
            r_slopes_norm[act_ix] = self.r_slopes[act_ix] / wt_r_avgs[act_ix]
        return r_slopes_norm

    def wt_release_avg_std(self):
        wt_r_avgs = np.mean(self.wt_r_slopes, axis=1)
        wt_r_stds = np.std(self.wt_r_slopes, axis=1, ddof=1)
        return (wt_r_avgs, wt_r_stds)

    def nbd_avg_std(self):
        n_slope_avgs = np.mean(self.n_slopes, axis=2)
        n_slope_stds = np.std(self.n_slopes, axis=2, ddof=1)
        return (n_slope_avgs, n_slope_stds)

    def nbd_norm_slopes(self):
        return (self.n_slopes / self.n_maxes) * 100

    def nbd_norm_avg_std(self):
        n_norm_slopes = self.nbd_norm_slopes()
        n_norm_slope_avgs = np.mean(n_norm_slopes, axis=2)
        n_norm_slope_stds = np.std(n_norm_slopes, axis=2, ddof=1)
        return (n_norm_slope_avgs, n_norm_slope_stds)

def calc_initial_rate_samples(df, nbd_sites, timepoint_ix=4,
                              normalize_nbd=False):
    """Get initial NBD/Tb release values from the data at the given timepoint.

    Returns an instance of InitialRateSamples containing the data at the
    given timepoint for all of the mutants/replicates.

    If normalize_nbd is set, then the raw NBD fluorescence data is converted to
    F/F0 values.
    """
    replicates = range(1, 4)
    activators = ['Bid', 'Bim']
    nbd_sites_filt = [s for s in nbd_sites if s != 'WT']
    # Lists for storing all of the various fits, slopes
    r_slopes = np.zeros((len(activators), len(nbd_sites_filt), len(replicates)))
    n_slopes = np.zeros((len(activators), len(nbd_sites_filt), len(replicates)))
    n_maxes = np.zeros((len(activators), len(nbd_sites_filt), len(replicates)))
    wt_r_slopes = np.zeros((len(activators), len(replicates)))
    # For both Bid and Bim...
    for act_ix, activator in enumerate(activators):
        # Get release data for WT as a reference
        if 'WT' in nbd_sites:
            for rep_index, rep_num in enumerate(replicates):
                ry = df[(activator, 'Release', 'WT', rep_num, 'VALUE')].values
                wt_r_slopes[act_ix, rep_index] = ry[timepoint_ix]
        # If we don't have an entry for WT, put a NaN into the matrix
        else:
            for rep_index, rep_num in enumerate(replicates):
                wt_r_slopes[act_ix, rep_index] = ry[timepoint_ix] = np.nan
        # Iterate over all of the mutants
        for nbd_index, nbd_site in enumerate(nbd_sites_filt):
            # Iterate over the replicates for this mutant. Note that rep_num is
            # the 1-indexed number of the replicate (1, 2, 3) whereas the index
            # is the 0-based index into the array for the replicates (0, 1, 2)
            for rep_index, rep_num in enumerate(replicates):
                # Get the release data
                rt = df[(activator, 'Release',
                         nbd_site, rep_num, 'TIME')].values
                ry = df[(activator, 'Release',
                         nbd_site, rep_num, 'VALUE')].values
                time = rt[timepoint_ix]
                # Store the y value at the given timepoint
                r_slopes[act_ix, nbd_index, rep_index] = ry[timepoint_ix]

                # Get the NBD data
                nt = df[(activator, 'NBD', nbd_site, rep_num, 'TIME')].values
                ny = df[(activator, 'NBD', nbd_site, rep_num, 'VALUE')].values
                # Normalize the NBD data
                if normalize_nbd:
                    ny = ny / ny[0]
                    f0 = 1.0
                else:
                    f0 = ny[0]
                # Maximum NBD F/F0
                # Because NBD-47-Bax actually decreases in fluorescence, we
                # use the change relative to the minimum.
                if nbd_site == '47':
                    n_slopes[act_ix, nbd_index, rep_index] = \
                                                        f0 - ny[timepoint_ix]
                    n_maxes[act_ix, nbd_index, rep_index] = f0 - np.min(ny)
                else:
                    n_slopes[act_ix, nbd_index, rep_index] = \
                                                        ny[timepoint_ix] - f0
                    n_maxes[act_ix, nbd_index, rep_index] = np.max(ny) - f0
    # Save the results in the initial rate structure
    irs = InitialRateSamples(r_slopes, n_slopes, n_maxes, wt_r_slopes, time)
    return irs

def plot_initial_rate_samples(df, nbd_sites, timepoint_ix=4,
                              file_basename=None, normalized_to_wt=True):
    """Plot characteristics of initial rates.
    """
    set_fig_params_for_publication()
    # Important lists
    replicates = range(1, 4)
    activators = ['Bid', 'Bim']
    nbd_sites_filt = [s for s in nbd_sites if s != 'WT']
    # Get the initial rate data (copy over to local variables for legacy
    # reasons)
    irs = calc_initial_rate_samples(df, nbd_sites, timepoint_ix,
                                    normalize_nbd=True)
    if normalized_to_wt:
        (r_slope_avgs, r_slope_stds) = irs.release_avg_std_wt_normalized()
        r_slopes = irs.release_slopes_wt_normalized()
    else:
        (r_slope_avgs, r_slope_stds) = irs.release_avg_std()
        r_slopes = irs.r_slopes

    (n_slope_avgs, n_slope_stds) = irs.nbd_avg_std()
    n_norm_slopes = irs.nbd_norm_slopes()
    (n_norm_slope_avgs, n_norm_slope_stds) = irs.nbd_norm_avg_std()
    time = irs.time
    # Bar plots of initial points ------------------------------------
    fig_names = {'Release': '%s_init_release_bar' % file_basename,
                 'NBD': '%s_init_nbd_bar' % file_basename}
    if normalized_to_wt:
        ylabels = {'Release': 'Dye release at %d sec\n(fold-change over WT)' %
                   time,
                   'NBD': '\\%% of Max NBD F/$F_0$ at %d sec' % time}
    else:
        ylabels = {'Release': '\\%% Dye release at %d sec' % time,
                   'NBD': '\\%% of Max NBD F/$F_0$ at %d sec' % time}
    # Make both Release and NBD bar plots
    for dtype in ['Release', 'NBD']:
        fig_name = fig_names[dtype]
        # Get the width of the figure based on the number of sites
        (fig_width, rel_left, rel_right) = \
                            calc_barplot_width(len(nbd_sites_filt))
        plt.figure(fig_name, figsize=(fig_width, 1.5), dpi=300)
        plt.ylabel(ylabels[dtype], fontsize=fontsize, multialignment='center')
        bar_colors = {'Bid': 'gray', 'Bim': 'black'}
        # Iterate over the activators
        for act_ix, activator in enumerate(activators):
            for nbd_index, nbd_site in enumerate(nbd_sites_filt):
                # Get the appropriate data
                if dtype == 'Release':
                    yval = r_slope_avgs[act_ix, nbd_index]
                    yerr = r_slope_stds[act_ix, nbd_index]
                elif dtype == 'NBD':
                    yval = n_norm_slope_avgs[act_ix, nbd_index]
                    yerr = n_norm_slope_stds[act_ix, nbd_index]
                else:
                    raise ValueError('Unknown datatype: %s' % dtype)
                # Plot the data
                plt.bar(range(nbd_index*3 + act_ix, (nbd_index*3) + 1 + act_ix),
                        yval, yerr=yerr, width=1, color=bar_colors[activator],
                        linewidth=0, ecolor='k', capsize=1.5)
        x_lbound = -1
        x_ubound = len(nbd_sites_filt) * 3
        # Add line showing wild type
        plt.hlines(0, x_lbound, x_ubound, linewidth=0.5)
        plt.subplots_adjust(left=rel_left, bottom=0.10, right=rel_right,
                            top=0.94)
        # Set bounds
        ax = plt.gca()
        ax.set_xlim([x_lbound, x_ubound])
        # Format the plot
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_tick_params(which='both', labelsize=fontsize, pad=2,
                length=2, width=0.5)
        ax.yaxis.set_tick_params(which='both', direction='out',
                        labelsize=fontsize, pad=0, length=2, width=0.5)
        ax.xaxis.labelpad = 2
        ax.yaxis.labelpad = 2
        ax.set_xticks(np.arange(1, 1 + len(nbd_sites_filt) * 3, 3))
        ax.set_xticklabels(nbd_sites_filt)
        # Create legend with rectangles to match colors of bars
        bid_patch = mpatches.Patch(color=bar_colors['Bid'], label='cBid')
        bim_patch = mpatches.Patch(color=bar_colors['Bim'], label='Bim')
        leg = plt.legend(handles=[bid_patch, bim_patch],
                bbox_to_anchor=(1.02, 0.5), loc='center left', borderaxespad=0.,
                prop={'size': fontsize}, handlelength=1)
        leg.draw_frame(False)
        # Output the file, if desired
        if file_basename:
            plt.savefig('%s.pdf' % fig_name)
            plt.savefig('%s.png' % fig_name)

    # 2D scatter plot of NBD/Tb rates, no normalization ----------------------
    for act_ix, activator in enumerate(activators):
        fig_name = "%s_init_scatter_no_norm_%s" % (file_basename, activator)
        plt.figure(fig_name, figsize=(2, 2), dpi=300)
        plt.errorbar(r_slope_avgs[act_ix], n_slope_avgs[act_ix],
                     xerr=r_slope_stds[act_ix], yerr=n_slope_stds[act_ix],
                     marker='o', color='k', markersize=3, linestyle='',
                     capsize=1.5)
        plt.xlabel('Dye release at %d sec' % time)
        plt.ylabel('delta F/$F_0$ at %d sec' % time)
        x_offset = 0.4
        y_offset = 0.04
        for nbd_index, nbd_site in enumerate(nbd_sites_filt):
            plt.text(r_slope_avgs[act_ix, nbd_index] + x_offset,
                     n_slope_avgs[act_ix, nbd_index] + y_offset, nbd_site,
                     fontsize=fontsize)
        ax = plt.gca()
        ybounds = ax.get_ylim()
        #plt.vlines(wt_avg, ybounds[0], ybounds[1], linestyle='--')
        format_axis(ax)
        plt.subplots_adjust(left=0.17, bottom=0.15)
        #if file_basename:
        #    plt.savefig('%s.pdf' % fig_name)

    # 2D scatter plot of NBD/Tb rates, normalized ----------------------------
    # Create lookup dict for sites
    site_region = dict(
         map(lambda s: (s, 'nterm'), ['3', '5', '15', '36', '40', '47']) +
         map(lambda s: (s, 'bh3'),  ['54', '62', '68', '79']) +
         map(lambda s: (s, 'h56'), ['120', '122', '126', '138', '151']) +
         map(lambda s: (s, 'h9'), ['175', '179', '184', '188']))
    # Color code the points based on the region of the protein
    color_dict = collections.OrderedDict([
                            ('nterm', 'purple'), ('bh3', 'red'),
                            ('h56', 'green'), ('h9', 'blue')])
    # Make two separate plots, one for each activator
    for act_ix, activator in enumerate(activators):
        fig_name = "%s_init_scatter_norm_%s" % (file_basename, activator)
        plt.figure(fig_name, figsize=(2, 2), dpi=300)
        plt.xlabel('Dye release at %d sec\n(fold-change over WT)' % time)
        plt.ylabel('\\%% of Max NBD F/$F_0$ at %d sec' % time)
        # Iterate over the sites, plotting one point at a time
        for nbd_index, nbd_site in enumerate(nbd_sites_filt):
            # Get color for site based on region
            pt_color = color_dict[site_region[nbd_site]]
            # Plot a point on the scatterplot
            plt.errorbar(r_slope_avgs[act_ix, nbd_index],
                         n_norm_slope_avgs[act_ix, nbd_index],
                         xerr=r_slope_stds[act_ix, nbd_index],
                         yerr=n_norm_slope_stds[act_ix, nbd_index],
                         marker='', color=pt_color, linestyle='',
                         markersize=2, capsize=1)
        #plt.vlines(wt_avg, np.min(n_slope_avgs), np.max(n_slope_avgs),
        #           linestyle='--')

        # Filter out the BH3 region data for linear fit
        x_pts = []
        y_pts = []
        for nbd_index, nbd_site in enumerate(nbd_sites_filt):
            # Skip residues in the BH3 region
            if site_region[nbd_site] == 'bh3':
                continue
            # If not in BH3 region, add to list of points
            for rep_ix, rep_num in enumerate(replicates):
                x_pts.append(r_slopes[act_ix, nbd_index, rep_ix])
                y_pts.append(n_norm_slopes[act_ix, nbd_index, rep_ix])
        # Turn the lists into numpy arrays
        x_pts = np.array(x_pts)
        y_pts = np.array(y_pts)
        # Fit the filtered data to a line, allowing free intercept
        linfit = scipy.stats.linregress(x_pts, y_pts)
        plt.plot(x_pts, x_pts * linfit[0] + linfit[1], color='k', zorder=1)
        print linfit

        # Format the plot
        ax = plt.gca()
        format_axis(ax)
        plt.subplots_adjust(left=0.22, bottom=0.19, right=0.95, top=0.95)
        ybounds = ax.get_ylim()
        ax.set_ylim(min(0, ybounds[0]), ybounds[1])
        # Create legend with lines to match colors of points
        n_patch = mlines.Line2D([], [], color=color_dict['nterm'],
                                label='N-term')
        bh3_patch = mlines.Line2D([], [], color=color_dict['bh3'],
                                  label='BH3')
        a56_patch = mlines.Line2D([], [], color=color_dict['h56'],
                                  label='$\\alpha$5-6')
        a9_patch = mlines.Line2D([], [], color=color_dict['h9'],
                                  label='$\\alpha 9$')
        leg = plt.legend(handles=[n_patch, bh3_patch, a56_patch, a9_patch],
                 loc='lower right', borderaxespad=0.,
                 prop={'size': fontsize}, handlelength=2,
                 handletextpad=0.2, labelspacing=0.5)
        leg.draw_frame(False)
        # Draw residue number text for BH3 residues only
        for nbd_index, nbd_site in enumerate(nbd_sites_filt):
            (x_lbound, x_ubound) = plt.xlim()
            x_range = x_ubound - x_lbound
            x_offset = x_range * -0.07
            y_offset = 1.15
            # Select residues
            if site_region[nbd_site] == 'bh3':
                # Get color for site text based on region
                pt_color = color_dict[site_region[nbd_site]]
                # Draw text
                plt.text(r_slope_avgs[act_ix, nbd_index] + x_offset,
                         n_norm_slope_avgs[act_ix, nbd_index] + y_offset,
                         nbd_site, color=pt_color, fontsize=fontsize)
        # Save the fig, if desired
        if file_basename:
            plt.savefig('%s.pdf' % fig_name)
            plt.savefig('%s.png' % fig_name)

    # Plot of all replicates
    for act_ix, activator in enumerate(activators):
        plt.figure(activator)
        for nbd_index, nbd_site in enumerate(nbd_sites_filt):
            for rep_ix, rep_num in enumerate(replicates):
                pt_color = color_dict[site_region[nbd_site]]
                # Plot point
                plt.plot(r_slopes[act_ix, nbd_index, rep_ix],
                         n_norm_slopes[act_ix, nbd_index, rep_ix],
                         marker='o', linestyle='', color=pt_color)
    return

def plot_rank_changes(means1, sds1, means2, sds2, nbd_sites):
    ubounds1 = means1 + sds1
    lbounds1 = means1 - sds1
    ubounds2 = means2 + sds2
    lbounds2 = means2 - sds2
    ubounds = np.vstack([ubounds1, ubounds2]).T
    lbounds = np.vstack([lbounds1, lbounds2]).T
    means = np.vstack([means1, means2]).T
    plt.figure()
    num_colors = ubounds.shape[0]
    colors = color_list = plt.cm.Set3(np.linspace(0, 1, num_colors))
    for i in range(ubounds.shape[0]):
        if nbd_sites[i] in ['175', '179', '184', '188']:
        #if nbd_sites[i] in ['3', '5', '15', '36', '40', '47']:
            color='k'
        else:
            color=colors[i]
        plt.fill_between([1, 2], ubounds[i], lbounds[i], color=color,
                         alpha=0.2)
        plt.plot([1, 2], means[i], color=color)
        plt.text(0.95, means[i, 0], nbd_sites[i], color='k',
                 fontsize=8)
        plt.text(2.05, means[i, 1], nbd_sites[i], color='k',
                 fontsize=8)
    plt.xlim(0.5, 2.5)

def plot_bid_vs_bim_release(df, nbd_sites, dtype='Release',
                                file_basename=None):
    replicates = range(1, 4)
    activators = ['Bid', 'Bim']
    # Get the length of the timecourses
    for nbd_index, nbd_site in enumerate(nbd_sites):
        color_ix = 0
        for act_ix, activator in enumerate(activators):
            # Initialization for WT
            wt_slice = df[activator][dtype]['WT']
            wt_numpts = wt_slice.shape[0]
            wt_y = np.zeros((wt_numpts, len(replicates)))
            # Initialization for mutant
            mut_slice = df[activator][dtype][nbd_site]
            mut_numpts = mut_slice.shape[0]
            mut_y = np.zeros((mut_numpts, len(replicates)))
            # Iterate over reps and get data
            for rep_ix, rep_num in enumerate(replicates):
                # Only get the time coordinates for the first rep
                if rep_ix == 0:
                    wt_time = wt_slice[rep_num]['TIME'].values
                    mut_time = mut_slice[rep_num]['TIME'].values
                wt_y[:, rep_ix] = wt_slice[rep_num]['VALUE'].values
                mut_y[:, rep_ix] = mut_slice[rep_num]['VALUE'].values
            # Now get the averages and SDs
            wt_avg = np.mean(wt_y, axis=1)
            wt_sd = np.std(wt_y, axis=1, ddof=1)
            wt_ubound = wt_avg + wt_sd
            wt_lbound = wt_avg - wt_sd
            mut_avg = np.mean(mut_y, axis=1)
            mut_sd = np.std(mut_y, axis=1, ddof=1)
            mut_ubound = mut_avg + mut_sd
            mut_lbound = mut_avg - mut_sd

            #num_colors = 4
            #colors = plt.cm.Set3(np.linspace(0, 1, num_colors))
            colors = ['r', 'm', 'b', 'g']
            fig_name = 'bid_bim_tc_comp_%s' % nbd_site
            # Plot the ratio
            plt.figure(fig_name, figsize=(10, 10))
            plt.subplot(1, 2, 1)
            (ratio_avg, ratio_sd) = \
                            calc_ratio_mean_sd(mut_avg, mut_sd, wt_avg, wt_sd)
            plt.plot(wt_time, ratio_avg, color=colors[color_ix],
                     label=activator)
            plt.fill_between(wt_time, ratio_avg + ratio_sd,
                             ratio_avg - ratio_sd, alpha=0.5,
                             color=colors[color_ix])
            plt.legend(loc='upper right', fontsize=10)
            plt.ylim(0, 5)

            # Plot the raw timecourses for WT and mutant
            # Plot the mutant
            plt.subplot(1, 2, 2)
            plt.plot(wt_time, mut_avg, color=colors[color_ix],
                     label='%s, NBD-%sC-Bax' % (activator, nbd_site))
            plt.fill_between(wt_time, mut_ubound, mut_lbound,
                             color=colors[color_ix], alpha=0.2)
            color_ix += 1
            # Plot the WT
            plt.plot(wt_time, wt_avg, color=colors[color_ix],
                     label='%s, WT Bax' % activator)
            plt.fill_between(wt_time, wt_ubound, wt_lbound,
                             color=colors[color_ix], alpha=0.3)
            plt.legend(loc='lower right', fontsize=10)
            color_ix += 1
            if file_basename:
                plt.savefig('%s_%s.pdf' % (file_basename, fig_name))
                plt.savefig('%s_%s.png' % (file_basename, fig_name))

def calc_release_peaks(df, nbd_sites, activators=None, replicates=None,
                       window=1, csv_filename=None):
    """Measure the lag phase of the release data.

    Takes the derivative of the release data for the given set of
    activators, NBD sites, and replicates, and gets the time associated
    with the peak of the derivative.

    Returns
    -------
    Dictionary containing keys of the form (activator, nbd_site, replicate)
    mapped to the times of the maximum rate of release.
    """

    # Set some defaults
    if activators is None:
        activators = ['Bid', 'Bim']
    if replicates is None:
        replicates = range(1, 4)
    peak_dict = collections.OrderedDict()
    # Initialize the filter
    b, a = scipy.signal.butter(1, 0.2)
    for activator, nbd_site, rep_index in \
            itertools.product(activators, nbd_sites, replicates):
        key = (activator, nbd_site, rep_index)
        # Get the data
        rt = df[(activator, 'Release', nbd_site,
                 rep_index, 'TIME')].values
        ry = df[(activator, 'Release', nbd_site,
                 rep_index, 'VALUE')].values
        # Apply the filter

        # Filter the timecourse
        r_filt = scipy.signal.filtfilt(b, a, ry)
        r_avg = moving_average(r_filt, n=window)
        # Take the derivative
        r_diff = np.diff(r_avg)
        # When we take the derivative, the size of the array shrinks by
        # one, because we are calculating the differences between neighboring
        # elements. So if the data is [0, 1, 3, 4, 5], the derivatives will
        # be [1, 2, 1, 1], and the maximum of the derivative array will be at
        # index 1, which corresponds to the difference of two and the entry of
        # three in the original array. If we adopt the convention that the
        # index to use for the maximum slope is the latter of the two values
        # used in calculating the difference, this means we need to add one to
        # the index associated with the maximum value of the diff array.
        r_max_tpt = np.argmax(r_diff) + 1
        peak_dict[key] = rt[r_max_tpt]

    if csv_filename:
        with open(csv_filename, 'w') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=',')
            for key, value in peak_dict.iteritems():
                line = list(key)
                line.append(value)
                csv_writer.writerow(line)
    return peak_dict

def plot_example_derivatives(df, activator, nbd_site, rep_index, window=1,
                             normalize_nbd=False, plot_filename=None):
    set_fig_params_for_publication()
    # Create an order 3 lowpass butterworth filter.
    b, a = scipy.signal.butter(1, 0.2)

    # RELEASE DERIVATIVE
    rt = df[(activator, 'Release', nbd_site, rep_index, 'TIME')].values
    ry = df[(activator, 'Release', nbd_site, rep_index, 'VALUE')].values
    # Filter the timecourse
    r_filt = scipy.signal.filtfilt(b, a, ry)
    r_avg = moving_average(r_filt, n=window)
    # Take the derivative
    r_diff = np.diff(r_avg)

    # NBD DERIVATIVE
    nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
    ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values
    # Normalize NBD to F/F0
    if normalize_nbd:
        ny = ny / float(ny[0])
    # Filter
    n_filt = scipy.signal.filtfilt(b, a, ny)
    n_avg = moving_average(n_filt, n=window)
    # Take derivative
    n_diff = np.diff(n_avg)

    # PLOT
    fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
    ax = fig.gca()
    n_diff_norm = n_diff / np.max(np.abs(n_diff))
    r_diff_norm = r_diff / np.max(np.abs(r_diff))
    ax.plot(nt[1+window-1:], n_diff_norm, label=r'$\frac{d}{dt}$ NBD')
    ax.plot(rt[1+window-1:], r_diff_norm, color='r',
            label=r'$\frac{d}{dt}$ Tb release')
    ax.set_xlabel('Time (sec)')
    ax.set_ylabel(r'\% Max Rate')
    #ax.set_title('%s, NBD-%s-Bax normalized derivative' % (activator, nbd_site))
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlim(0, 2000)
    plt.subplots_adjust(left=0.22, bottom=0.19)
    plt.legend(loc='upper right', fontsize=fontsize, frameon=False)
    format_axis(ax)

    if plot_filename:
        plt.savefig('%s.pdf' % plot_filename)
        plt.savefig('%s.png' % plot_filename, dpi=300)

def plot_derivatives(df, nbd_sites, normalize_nbd=False):
    replicates = range(1, 4)
    num_pts = 4
    window = 1 # For moving average
    activators = ['Bid', 'Bim']
    # Create an order 3 lowpass butterworth filter.
    b, a = scipy.signal.butter(1, 0.2)
    for nbd_index, nbd_site in enumerate(nbd_sites):
        for activator in activators:
            # We store the list of timepoints where the release derivative
            # reaches its peak so that we can plot lines for all three with
            # the same upper and lower y-coordinates.
            peak_pts = []
            for rep_index in replicates:
                rt = df[(activator, 'Release', nbd_site,
                         rep_index, 'TIME')].values
                ry = df[(activator, 'Release', nbd_site,
                         rep_index, 'VALUE')].values
                # Filter the timecourse
                r_filt = scipy.signal.filtfilt(b, a, ry)
                r_avg = moving_average(r_filt, n=window)
                # Take the derivative
                r_diff = np.diff(r_avg)
                # See comment in calc_release_peaks, above
                r_max_tpt = np.argmax(r_diff) + 1
                peak_pts.append(rt[r_max_tpt])
                # Calculate max NBD slope, but not for WT
                if nbd_site != 'WT':
                    nt = df[(activator, 'NBD', nbd_site,
                             rep_index, 'TIME')].values
                    ny = df[(activator, 'NBD', nbd_site,
                             rep_index, 'VALUE')].values
                    # Normalize NBD to F/F0
                    if normalize_nbd:
                        ny = ny / float(ny[0])
                    # Filter
                    n_filt = scipy.signal.filtfilt(b, a, ny)
                    n_avg = moving_average(n_filt, n=window)
                    # Take derivative
                    n_diff = np.diff(n_avg)

                # Terbium subplot
                plt.figure('%s, NBD-%s-Bax derivative' % (activator, nbd_site),
                           figsize=(12, 5))
                plt.subplot(1, 2, 1)
                plt.plot(rt[1+window-1:], r_diff, color=rep_colors[rep_index],
                         label='%s Rep %d' % (activator, rep_index))
                plt.ylabel('dRel/dt (% rel $sec^{-1}$)')
                plt.title('%s, NBD-%s-Bax, Tb derivative' %
                          (activator, nbd_site))
                plt.legend(loc='upper right')

                if nbd_site != 'WT':
                    # NBD subplot
                    plt.subplot(1, 2, 2)
                    plt.plot(nt[1+window-1:], n_diff,
                             color=rep_colors[rep_index],
                             label='%s Rep %d' % (activator, rep_index))
                    plt.xlabel('Time (sec)')
                    plt.ylabel('dNBD/dt ($F/F_0\ sec^{-1}$)')
                    plt.title('%s, NBD-%s-Bax, NBD derivative' %
                              (activator, nbd_site))
                    plt.legend(loc='upper right')

                    # Plot normalized derivatives
                    plt.figure('%s, NBD-%s-Bax normalized derivative' %
                               (activator, nbd_site))
                    n_diff_norm = n_diff / np.max(np.abs(n_diff))
                    r_diff_norm = r_diff / np.max(np.abs(r_diff))
                    plt.plot(rt[1+window-1:], r_diff_norm,
                             color=rep_colors[rep_index],
                             linestyle=line_styles[2],
                             label='%s Rep %d' % (activator, rep_index))
                    plt.plot(nt[1+window-1:], n_diff_norm,
                             color=rep_colors[rep_index],
                             linestyle=line_styles[3])
                    plt.xlabel('Time (sec)')
                    plt.ylabel('% max rate')
                    plt.title('%s, NBD-%s-Bax normalized derivative' %
                              (activator, nbd_site))
                    plt.legend(loc='upper right')

                # Add vertical lines to the normalized derivative plot
                plt.figure('%s, NBD-%s-Bax normalized derivative' %
                           (activator, nbd_site))
                ymin, ymax = plt.ylim()
                plt.vlines(peak_pts, ymin, ymax)

                # Call tight_layout for the Tb/NBD 2-panel figure
                plt.figure('%s, NBD-%s-Bax derivative' % (activator, nbd_site))
                plt.tight_layout()


def welch_t_test(means1, sds1, means2, sds2):
    n1 = 3
    n2 = 3
    t_numer = means1 - means2
    sq_sum = ((sds1**2)/n1) + ((sds2**2)/n2)
    t_denom = np.sqrt(sq_sum)
    t = t_numer / t_denom
    print t

    v_numer = sq_sum ** 2
    v1 = n1 - 1.0
    v2 = n2 - 1.0
    v_denom = ((sds1 ** 4) / ((n1**2) * v1)) / ((sds2 ** 4) / ((n2**2) * v2))
    v = v_numer / v_denom
    print v

    p_val = scipy.stats.t.sf(t, v)
    print p_val

def student_t_test(means1, sds1, means2, sds2, n):
    n = float(n)
    t_numer = np.abs(means1 - means2)
    sx1x2 = np.sqrt(0.5*(sds1**2 + sds2**2))
    t_denom = sx1x2 * np.sqrt(2/n)
    t = t_numer / t_denom
    print t

    p_val = scipy.stats.t.sf(t, n - 1)
    print p_val * 2
    return p_val * 2

# -- Model fits --
params_dict = {'c1_to_c2_k': 1e-4, 'c1_scaling': 2,
               'c0_to_c1_k': 2e-3}

def plot_2conf_fits(df, nbd_sites, activator, normalize_nbd=False):
    replicates = range(1, 4)
    fit_results = []
    # Filter out the WT residue, if present in the list
    nbd_sites_filt = [s for s in nbd_sites if s != 'WT']

    for nbd_index, nbd_site in enumerate(nbd_sites_filt):
        k1_means = []
        k1_sds = []

        for rep_index in replicates:
            nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
            ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values
            # Normalize NBD to F/F0
            if normalize_nbd:
                ny = ny / ny[0]

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

        plt.figure('%s, NBD-%s-Bax Fits' % (activator, nbd_site))
        plt.tight_layout()

        plt.figure("Fitted k1")
        plt.bar(range(nbd_index*4, (nbd_index*4) + 3), k1_means,
                yerr=k1_sds, width=1, color='r', ecolor='k')

    num_sites = len(nbd_sites_filt)
    plt.figure("Fitted k1")
    plt.ylabel('Fitted k1 ($sec^{-1}$)')
    ax = plt.gca()
    ax.set_xticks(np.arange(1.5, 1.5 + num_sites * 4, 4))
    ax.set_xticklabels(nbd_sites_filt)

    return fit_results

def plot_3conf_fits(df, nbd_sites, activator, normalize_nbd=False):
    #nbd_sites = ['15', '54', '62', '68', '79', '126', '138', '175']
    # Filter out the WT residue, if present in the list
    nbd_sites_filt = [s for s in nbd_sites if s != 'WT']
    replicates = range(1, 4)
    fit_results = []
    for nbd_index, nbd_site in enumerate(nbd_sites_filt):
        k1_means = []
        k2_means = []
        k1_sds = []
        k2_sds = []

        for rep_index in replicates:
            rt = df[(activator, 'Release', nbd_site, rep_index, 'TIME')].values
            ry = df[(activator, 'Release', nbd_site, rep_index, 'VALUE')].values
            nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
            ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values
            # Normalize NBD to F/F0
            if normalize_nbd:
                ny = ny / ny[0]
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

    num_sites = len(nbd_sites_filt)
    plt.figure("Fitted k1/k2")
    plt.ylabel('Fitted k1/k2 ($sec^{-1}$)')
    ax = plt.gca()
    ax.set_xticks(np.arange(3, 3 + num_sites * 7, 7))
    ax.set_xticklabels(nbd_sites_filt)

    return fit_results

def plot_nbd_error_estimates(df, nbd_sites, dtype='NBD', activators=None,
                             replicates=None, last_n_pts=50, fit_type='cubic',
                             plot=True, normalize_nbd=False):
    # Filter out the WT residue, if present in the list
    nbd_sites_filt = [s for s in nbd_sites if s != 'WT']
    # Set some defaults
    if activators is None:
        activators = ['Bid', 'Bim']
    if replicates is None:
        replicates = (1, 2, 3)

    for nbd_index, nbd_site in enumerate(nbd_sites_filt):
        for act_ix, activator in enumerate(activators):
            residuals_reps = []
            for rep_index, rep_num in enumerate(replicates):
                ny = df[activator, dtype, nbd_site, rep_num, 'VALUE'].values
                # Normalize NBD to F/F0
                if dtype == 'NBD' and normalize_nbd:
                    ny = ny / ny[0]
                # Get the error value and the fig (if not figure desired, then
                # fig will be None
                (residuals, fig) = cev.calc_err_var(ny, last_n_pts=last_n_pts,
                                           fit_type=fit_type, plot=plot)
                residuals_reps.append(residuals)
                # If we're plotting, we should have a figure object, and if
                # not, then not
                assert ((plot is None) == (fig is None))
                # Add title info to plot
                nbd_err = np.std(residuals, ddof=1)
                if plot:
                    fig.subplots_adjust(top=0.86)
                    fig.text(0.5, 0.95,
                             'Est. Error for %s, %s, %sC-Bax, rep %d: %f' %
                             (activator, dtype, nbd_site, rep_num, nbd_err),
                             verticalalignment='top',
                             horizontalalignment='center')

            # If there is more than one replicates, combine the three residuals
            # replicates into one big set of residuals
            if len(replicates) == 1:
                continue
            pooled_residuals = np.array(residuals_reps).flatten()
            pooled_err = np.std(pooled_residuals, ddof=1)
            if plot:
                fig = plt.figure(figsize=(8, 5))
                # Plot of fit and distribution of residuals
                plt.subplot(1, 2, 1)
                plt.hist(pooled_residuals)
                plt.title('Histogram of pooled residuals')
                plt.subplot(1, 2, 2)
                scipy.stats.probplot(pooled_residuals, dist='norm', plot=plt)
                plt.title('Quantile-quantile plot vs. normal')
                plt.tight_layout()
                plt.subplots_adjust(top=0.86)
                fig.text(0.5, 0.95,
                         'Est. Error of %s for %s, %sC-Bax: %f' %
                         (activator, dtype, nbd_site, pooled_err),
                         verticalalignment='top',
                         horizontalalignment='center')

if __name__ == '__main__':
    from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
    #from tbidbaxlipo.data.parse_bid_bim_fret_nbd_release import df, nbd_residues
    plt.ion()
    #import pickle
    #with open('data2.pck', 'w') as f:
        #pickle.dump((df, nbd_residues), f)
    #with open('data2.pck') as f:
    #    (df, nbd_residues) = pickle.load(f)
    #plot_release_endpoints(df, nbd_residues, normalized_to_wt=True,
    #                       last_n_pts=3, file_basename=None)
    #plot_initial_rate_samples(df, nbd_residues, timepoint_ix=20,
    #                          file_basename=None, normalized_to_wt=True)
    plot_example_derivatives(df, 'Bid', '15', 1)
    #plot_derivatives(df, ['WT'])
    #plot_nbd_error_estimates(df, ['68'], last_n_pts=80, fit_type='cubic')
    #plot_release_endpoints(df, nbd_residues, normalized_to_wt=True)
    #plot_bid_vs_bim_release(df, nbd_residues)
    #plot_nbd_error_estimates(df, ['54', '62', '68'])

    #irs = calc_initial_rate_samples(df, nbd_residues, timepoint_ix=15)
    #(r_slope_avgs, r_slope_stds) = irs.release_avg_std()
    #(n_slope_avgs, n_slope_stds) = irs.nbd_avg_std()
    #nbd_sites_filt = [s for s in nbd_residues if s != 'WT']
    #plot_rank_changes(r_slope_avgs[0], r_slope_stds[0],
    #                  r_slope_avgs[1], r_slope_stds[1], nbd_sites_filt)
    #plot_rank_changes(n_slope_avgs[0], n_slope_stds[0],
    #                  n_slope_avgs[1], n_slope_stds[1], nbd_sites_filt)
    #(ra, rs) = irs.release_avg_std_wt_normalized()
    #p_vals = student_t_test(ra[0], rs[0], ra[1], rs[1], 3)
    #print p_vals < 0.1
