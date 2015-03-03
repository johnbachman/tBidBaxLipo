import collections
from matplotlib import pyplot as plt
import numpy as np
import scipy.stats
from tbidbaxlipo.util import set_fig_params_for_publication, fontsize, \
                             format_axis
from tbidbaxlipo.util.error_propagation import calc_ratio_mean_sd
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from tbidbaxlipo.util import fitting

line_colors = {'Bid': 'r', 'Bim': 'b'}
dtype_line_colors = {'Release':'r', 'NBD':'g', 'FRET':'b'}
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

def plot_all(df, nbd_residues, datatypes, file_basename=None):
    """Release should be first in the list of datatypes."""
    if datatypes[0] != 'Release':
        raise ValueError("Release should be first in the datatypes list.")
    # Some definitions
    activators = ['Bid', 'Bim']
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
        plt.figure(figsize=(fig_width, 5))
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
                for i in range(1, 4):
                    # Get the data
                    t = df[(activator, dtype, nbd_site, i, 'TIME')]
                    v = df[(activator, dtype, nbd_site, i, 'VALUE')]
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
        plt.tight_layout()
        if file_basename:
            plt.savefig('%s_%s.pdf' % (file_basename, nbd_site))

def plot_all_by_replicate(df, nbd_residues, datatypes, file_basename=None):
    """Release should be first in the datatypes list."""
    if datatypes[0] != 'Release':
        raise ValueError("Release should be first in the datatypes list.")
    # Some definitions
    activators = ['Bid', 'Bim']
    # Every mutant gets its own plot
    for nbd_index, nbd_site in enumerate(nbd_residues):
        for activator in activators:
            plt.figure(figsize=(14, 5))
            # Each replicate gets its own subplot
            ax2_list = []
            nbd_max = -np.inf
            nbd_min = np.inf
            for rep in [1, 2, 3]:
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

def calc_barplot_width(num_sites):
    # I found that these numbers worked for a 4 inch wide figure with the
    # given relative left and right margins and 19 sites. This allows the
    # same proportions for endpoint plots with different numbers of sites.
    rmarg = 0.13
    lmarg = 0.11
    abs_left = 4 * lmarg
    abs_right = 4 * rmarg
    abs_middle = 4 * (1 - lmarg - rmarg)
    abs_per_site = abs_middle / 19.
    fig_width = abs_per_site * num_sites + abs_left + abs_right
    rel_left = abs_left / float(fig_width)
    rel_right = 1 - (abs_right / float(fig_width))
    return (fig_width, rel_left, rel_right)

def plot_nbd_endpoints(df, nbd_sites, last_n_pts=3, file_basename=None):
    replicates = range(1, 4)
    activators = ['Bid', 'Bim']
    # Filter out the WT residue from the list, if its there
    nbd_sites_no_wt = [s for s in nbd_sites if s != 'WT']
    # Matrix for storing the endpoints for all mutants, replicates
    n_endpts = np.zeros((len(nbd_sites_no_wt), len(replicates)))
    # Figure setup
    set_fig_params_for_publication()
    (fig_width, rel_left, rel_right) = \
                        calc_barplot_width(len(nbd_sites_no_wt))
    plt.figure(file_basename, figsize=(fig_width, 1.5), dpi=300)
    plt.ylabel('NBD F/$F_0$', fontsize=fontsize)
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
                n_endpts[nbd_index, rep_index] = np.mean(ny[-last_n_pts:])

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
    plt.hlines(range(20, 120, 20), x_lbound, x_ubound, color='0.85',
               zorder=1, linewidth=0.5)
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
    # Output the file, if desired
    if file_basename:
        plt.savefig('%s.pdf' % file_basename)

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

def calc_initial_rate_samples(df, nbd_sites, timepoint_ix=4):
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
                # Maximum NBD F/F0
                # Because NBD-47-Bax actually decreases in fluorescence, we
                # use the change relative to the minimum.
                if nbd_site == '47':
                    n_slopes[act_ix, nbd_index, rep_index] = \
                                                        1 - ny[timepoint_ix]
                    n_maxes[act_ix, nbd_index, rep_index] = 1 - np.min(ny)
                else:
                    n_slopes[act_ix, nbd_index, rep_index] = \
                                                        ny[timepoint_ix] - 1
                    n_maxes[act_ix, nbd_index, rep_index] = np.max(ny) - 1
    # Save the results in the initial rate structure
    irs = InitialRateSamples(r_slopes, n_slopes, n_maxes, wt_r_slopes, time)
    return irs

def plot_initial_rate_samples(df, nbd_sites, timepoint_ix=4,
                              file_basename=None, normalized_to_wt=True):
    set_fig_params_for_publication()
    # Important lists
    replicates = range(1, 4)
    activators = ['Bid', 'Bim']
    nbd_sites_filt = [s for s in nbd_sites if s != 'WT']
    # Get the initial rate data (copy over to local variables for legacy
    # reasons)
    irs = calc_initial_rate_samples(df, nbd_sites, timepoint_ix)
    if normalized_to_wt:
        (r_slope_avgs, r_slope_stds) = irs.release_avg_std_wt_normalized()
    else:
        (r_slope_avgs, r_slope_stds) = irs.release_avg_std()

    (n_slope_avgs, n_slope_stds) = irs.nbd_avg_std()
    n_norm_slopes = irs.nbd_norm_slopes()
    (n_norm_slope_avgs, n_norm_slope_stds) = irs.nbd_norm_avg_std()
    r_slopes = irs.r_slopes
    time = irs.time

    # Bar plots of initial points ------------------------------------
    fig_names = {'Release': '%s_init_release_bar' % file_basename,
                 'NBD': '%s_init_nbd_bar' % file_basename}
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

    # 2D scatter plot of NBD/Tb rates, no normalization ----------------------
    for act_ix, activator in enumerate(activators):
        fig_name = "%s_init_scatter_no_norm_%s" % (file_basename, activator)
        plt.figure(fig_name, figsize=(2, 2), dpi=300)
        plt.errorbar(r_slope_avgs[act_ix], n_slope_avgs[act_ix],
                     xerr=r_slope_stds[act_ix], yerr=n_slope_stds[act_ix],
                     marker='o', color='k', markersize=3, linestyle='',
                     capsize=1.5)
        plt.xlabel('\\%% Dye release at %d sec' % time)
        plt.ylabel('\\%% of Max NBD F/$F_0$ at %d sec' % time)
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
        #plt.subplots_adjust(left=rel_left, bottom=0.10, right=rel_right,
        #                    top=0.94)
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
        plt.xlabel('\\%% Dye release at %d sec' % time)
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
        plt.subplots_adjust(left=0.22, bottom=0.15, right=0.95, top=0.95)
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

if __name__ == '__main__':
    from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
    plt.ion()
    plot_release_endpoints(df, nbd_residues, normalized_to_wt=True)
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
