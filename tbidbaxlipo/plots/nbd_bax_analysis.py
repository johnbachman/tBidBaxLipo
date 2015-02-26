from matplotlib import pyplot as plt

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

def plot_all(df, nbd_residues, file_basename=None):
    for nbd_index, nbd_site in enumerate(nbd_residues):
        plt.figure(figsize=(11, 5))
        # Make the release plot
        plt.subplot(1, 2, 1)
        activators = ['Bid', 'Bim']
        for activator in activators:
            for i in range(1, 4):
                t = df[(activator, 'Release', nbd_site, i, 'TIME')]
                v = df[(activator, 'Release', nbd_site, i, 'VALUE')]
                plt.plot(t, v, label='%s Rep %d' % (activator, i),
                        color=line_colors[activator],
                        linestyle=line_styles[i])

                plt.xlabel('Time (sec)')
                plt.ylabel('Pct. Release')
                plt.title('Release for NBD-%s-Bax' % nbd_site)
                plt.legend(loc='lower right')
        # There is no NBD curve for WT Bax, so skip the NBD
        # plot
        if nbd_site == 'WT':
            continue
        # Make the NBD plot
        plt.subplot(1, 2, 2)
        for activator in ['Bid', 'Bim']:
            for i in range(1, 4):
                t = df[(activator, 'NBD', nbd_site, i, 'TIME')]
                v = df[(activator, 'NBD', nbd_site, i, 'VALUE')]
                plt.plot(t, v, label='%s Rep %d' % (activator, i),
                        color=line_colors[activator],
                        linestyle=line_styles[i])
                plt.xlabel('Time (sec)')
                plt.ylabel('F/F0')
                plt.title('F/F0 for NBD-%s-Bax' % nbd_site)
                plt.legend(loc='lower right')
        plt.tight_layout()
        if file_basename:
            plt.savefig('%s_%s.pdf' % (file_basename, nbd_index))
            plt.savefig('%s_%s.png' % (file_basename, nbd_index))

def plot_endpoints(df, nbd_sites, last_n_pts=3, file_basename=None):
    replicates = range(1, 4)
    activators = ['Bid', 'Bim']
    # Filter out the WT residue from the list, if its there
    nbd_sites_no_wt = [s for s in nbd_sites if s != 'WT']
    # Matrix for storing the endpoints for all mutants, replicates
    r_endpts = np.zeros((len(nbd_sites_no_wt), len(replicates)))
    n_endpts = np.zeros((len(nbd_sites_no_wt), len(replicates)))
    # Figure setup
    set_fig_params_for_publication()
    plt.figure('release_endpt', figsize=(4, 1.5), dpi=300) # Release figure
    plt.ylabel(r'\% Dye Release' + '\n(normalized to WT)',
               fontsize=fontsize, multialignment='center')

    plt.figure('nbd_endpt', figsize=(4, 1.5), dpi=300) # NBD figure
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
        for nbd_index, nbd_site in enumerate(nbd_sites_no_wt):
            if nbd_site == 'WT':
                continue
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
            plt.figure('release_endpt')
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
            plt.figure('nbd_endpt')
            plt.bar(range(nbd_index*3 + act_ix, (nbd_index*3) + 1 + act_ix),
                    np.mean(n_endpts[nbd_index, :]),
                    width=1, color=bar_colors[activator], linewidth=0,
                    ecolor='k', capsize=1.5,
                    yerr=np.std(n_endpts[nbd_index, :], ddof=1))

    fig_names = ['release_endpt', 'nbd_endpt']
    for fig_name in fig_names:
        plt.figure(fig_name)
        plt.subplots_adjust(left=0.11, bottom=0.10, right=0.97, top=0.94)
        ax = plt.gca()
        ax.set_xlim([0, len(nbd_sites_no_wt) * 3])
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
        plt.savefig('%s_%s.pdf' % (file_basename, fig_name))


if __name__ == '__main__':
    from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
    plot_all(df, nbd_residues)
