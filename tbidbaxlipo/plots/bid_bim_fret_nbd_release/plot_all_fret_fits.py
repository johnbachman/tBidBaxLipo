from tbidbaxlipo.util import set_fig_params_for_publication, format_axis, \
                             fontsize
from tbidbaxlipo.plots.bid_bim_fret_nbd_release.preprocess_data import \
        nbd_residues, df_pre as df
import numpy as np
from matplotlib import pyplot as plt
import itertools
import cPickle
import sys
import os

set_fig_params_for_publication()

num_samples = 100
module_name = 'tbidbaxlipo.plots.bid_bim_fret_nbd_release'
#mcmc_path = os.path.join(os.path.dirname(sys.modules[module_name].__file__),
#                         'fret_mcmc')
mcmc_dir = sys.argv[1]
mcmc_path = os.path.join(os.path.dirname(__file__), mcmc_dir)

norm = True if mcmc_dir.endswith('norm') else False

plot_path = 'data2_fit_plots'
if norm:
    plot_path += '_norm'

activators = ['Bid', 'Bim']
#nbd_residues = ['3', '54', '126']

for res in nbd_residues:
    for act in activators:
        fig, axarr = plt.subplots(2, 3, sharex=True, figsize=(4.5, 4.5),
                                  dpi=300)
        for col_ix, rep in enumerate((1, 2, 3)):
            if norm:
                filename = os.path.join(mcmc_path,
                            'pt_data2_fret_norm_%s_NBD_%s_r%s_3confs.mcmc' %
                            (act, res, rep))
            else:
                filename = os.path.join(mcmc_path,
                            'pt_data2_fret_%s_NBD_%s_r%s_3confs.mcmc' %
                            (act, res, rep))
            with open(filename) as f:
                print("Loading %s" % filename)
                (gf, sampler) = cPickle.load(f)
            for row_ix, obs in enumerate(['NBD', 'FRET']):
                # Do NBD plot
                ax1 = axarr[row_ix, col_ix]
                # Get data
                t = df[(act, obs, res, rep, 'TIME')].values
                v = df[(act, obs, res, rep, 'VALUE')].values
                ax1.plot(t, v, color='k')
                # Format axes
                act_name = 'cBid' if act == 'Bid' else act
                ax1.set_title('NBD-%sC-Bax Rep %d, DAC-%s' %
                              (res, rep, act_name), fontsize=fontsize)
                ax1.set_xlim(-50, 4000)
                ax1.set_xticks(np.linspace(0, 4000, 5))
                # Label the axes
                if row_ix == 1:
                    ax1.set_xlabel(r'Time (sec $\times 10^{-3}$)')
                    ax1.set_xticklabels([int(n) for n in np.linspace(0, 4, 5)])
                if obs == 'NBD':
                    ax1.set_ylabel('NBD F/$F_0$')
                else:
                    ax1.set_ylabel('FRET')

                num_params = sampler.flatchain.shape
                # Set obs_func
                gf.builder.set_nbd_fret_obs_func()
                # Sample over chain
                (nwalkers, nsteps) = sampler.chain.shape[1:3]
                walker_indices = np.random.randint(0, nwalkers,
                                                   size=num_samples)
                step_indices = np.random.randint(0, nsteps, size=num_samples)
                line_color = 'r' if obs == 'NBD' else 'g'
                for i in xrange(num_samples):
                    psamp = sampler.chain[0, walker_indices[i],
                                          step_indices[i], :]
                    # Set values of parameters
                    for p_ix, param in enumerate(gf.builder.estimate_params):
                        param.value = 10 ** psamp[p_ix]
                    ysim = gf.builder.obs_func(t)[obs]
                    ax1.plot(t, ysim, color=line_color, alpha=0.1) 

                """
                # Plot max a posterior fit
                maxp_flat_ix = np.argmax(sampler.lnprobability[0])
                maxp_ix = np.unravel_index(maxp_flat_ix,
                                           sampler.lnprobability[0].shape)
                pmax = sampler.chain[0, maxp_ix[0], maxp_ix[1]]
                # Set values of parameters
                for p_ix, param in enumerate(gf.builder.estimate_params):
                    param.value = 10 ** pmax[p_ix]
                ysim = gf.builder.obs_func(t)[obs]
                ax1.plot(t, ysim, color='m')
                """

                # Format axis
                format_axis(ax1)

        plt.subplots_adjust(left=0.14, bottom=0.14, hspace=0.35, wspace=0.4)
        if norm:
            filebase = '%s/pt_data2_fret_norm_%s_%s_3confs_fits' % \
                    (plot_path, act, res)
        else:
            filebase = '%s/pt_data2_fret_%s_%s_3confs_fits' % \
                    (plot_path, act, res)
        fig.savefig('%s.pdf' % filebase)
        fig.savefig('%s.png' % filebase)
        #fig.savefig('data3_example_nbd_fret_fits.pdf', dpi=300)
        #fig.savefig('data3_example_nbd_fret_fits.png', dpi=300)


