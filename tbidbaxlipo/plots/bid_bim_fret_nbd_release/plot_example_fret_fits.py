from tbidbaxlipo.util import set_fig_params_for_publication, format_axis, \
                             fontsize
from tbidbaxlipo.plots.bid_bim_fret_nbd_release.preprocess_data import \
        df_pre as df
import numpy as np
from matplotlib import pyplot as plt
import cPickle
import tbidbaxlipo.plots.bid_bim_fret_nbd_release
import sys
import os

set_fig_params_for_publication()

curves_to_plot = [('Bid', '54', 'NBD'),
                  ('Bid', '126', 'FRET')]
num_samples = 100
module_name = 'tbidbaxlipo.plots.bid_bim_fret_nbd_release'
mcmc_path = os.path.join(os.path.dirname(sys.modules[module_name].__file__),
                         'fret_mcmc')

fig, axarr = plt.subplots(2, 2, sharex=True, figsize=(3, 3), dpi=300)

for col_ix, curve_info in enumerate(curves_to_plot):
    act = curve_info[0]
    res = curve_info[1]
    for row_ix, obs in enumerate(['NBD', 'FRET']):
        # Do NBD plot
        ax1 = axarr[row_ix, col_ix]
        # Get data
        t = df[(act, obs, res, 1, 'TIME')].values
        v = df[(act, obs, res, 1, 'VALUE')].values
        ax1.plot(t, v, color='k')
        # Format axes
        act_name = 'cBid' if act == 'Bid' else act
        ax1.set_title('NBD-%sC-Bax, DAC-%s' % (res, act_name),
                      fontsize=fontsize)
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

        filename = os.path.join(mcmc_path,
                        'pt_data2_fret_%s_NBD_%s_r1_3confs.mcmc' % (act, res))
        with open(filename) as f:
            print("Loading %s" % filename)
            (gf, sampler) = cPickle.load(f)

        num_params = sampler.flatchain.shape
        # Set obs_func
        gf.builder.set_nbd_fret_obs_func()
        # Sample over chain
        (nwalkers, nsteps) = sampler.chain.shape[1:3]
        walker_indices = np.random.randint(0, nwalkers, size=num_samples)
        step_indices = np.random.randint(0, nsteps, size=num_samples)
        line_color = 'r' if obs == 'NBD' else 'g'
        for i in xrange(num_samples):
            psamp = sampler.chain[0, walker_indices[i], step_indices[i], :]
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
fig.savefig('data2_example_nbd_fret_fits.pdf', dpi=300)
fig.savefig('data2_example_nbd_fret_fits.png', dpi=300)


