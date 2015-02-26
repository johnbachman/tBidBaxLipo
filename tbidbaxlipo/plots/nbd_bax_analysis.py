from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import set_fig_params_for_publication, fontsize
from tbidbaxlipo.util.error_propagation import calc_ratio_mean_sd

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

def plot_all(df, nbd_residues, datatypes, file_basename=None):
    """Release should be first in the list of datatypes."""
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
            plt.savefig('%s_%s.pdf' % (file_basename, nbd_index))
            plt.savefig('%s_%s.png' % (file_basename, nbd_index))

def plot_nbd_endpoints(df, nbd_sites, last_n_pts=3, file_basename=None):
    replicates = range(1, 4)
    activators = ['Bid', 'Bim']
    # Filter out the WT residue from the list, if its there
    nbd_sites_no_wt = [s for s in nbd_sites if s != 'WT']
    # Matrix for storing the endpoints for all mutants, replicates
    n_endpts = np.zeros((len(nbd_sites_no_wt), len(replicates)))
    # Figure setup
    set_fig_params_for_publication()
    plt.figure(file_basename, figsize=(4, 1.5), dpi=300) # NBD figure
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
    plt.subplots_adjust(left=0.11, bottom=0.10, right=0.97, top=0.94)
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
    plt.figure(file_basename, figsize=(4, 1.5), dpi=300) # Release figure
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
    plt.subplots_adjust(left=0.11, bottom=0.10, right=0.97, top=0.94)
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
    if file_basename:
        plt.savefig('%s.pdf' % file_basename)
        plt.savefig('%s.png' % file_basename)


if __name__ == '__main__':
    from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
    plot_all(df, nbd_residues)
