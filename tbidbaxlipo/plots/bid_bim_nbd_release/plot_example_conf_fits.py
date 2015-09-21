from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis, \
                             fontsize
import sys
import cPickle

set_fig_params_for_publication()

if len(sys.argv) < 2:
    usage = 'Usage: python plot_example_2conf_fits.py num_confs'
    print usage
    sys.exit(1)
num_confs = int(sys.argv[1])

#fig = plt.figure(figsize=(3, 2), dpi=300)

curves_to_plot = [
        ('3', 1),
        #('5', 1),
        #('15', 1),
        ('47', 1),
        ('54', 2),
        #('62', 2),
        ('68', 1),
        ('120', 1),
        ('175', 1),
        ]
# Get data for NBD-54C-Bax
nsamples = 100
end_ix = -1
fig, axarr = plt.subplots(2, 3, sharex=True, figsize=(3, 2), dpi=300)
for plot_index, plot_tuple in enumerate(curves_to_plot):
    residue, rep = plot_tuple
    # Get the current axis
    (row_ix, col_ix) = np.unravel_index(plot_index, axarr.shape)
    ax = axarr[row_ix, col_ix]

    # Get the data from the dataframe
    t = df[('Bid', 'NBD', residue, rep, 'TIME')][:end_ix]
    v = df[('Bid', 'NBD', residue, rep, 'VALUE')][:end_ix]
    # Normalize the NBD values
    norm_v = v / float(v[0])
    # Plot the data
    ax.plot(t, norm_v, label=residue, color='r', linestyle='', markersize=2,
            marker='.')

    # Load the GlobalFit and MCMC chain
    mcmc_filename = 'mcmc/pt_data1_Bid_NBD_%s_r%s_%sconfs.mcmc' % \
                    (residue, rep, num_confs)
    with open(mcmc_filename) as f:
        (gf, sampler) = cPickle.load(f)
    # Get the randomized sample indices
    (nwalkers, nsteps) = sampler.chain.shape[1:3]
    walker_indices = np.random.randint(0, nwalkers, size=nsamples)
    step_indices = np.random.randint(0, nsteps, size=nsamples)
    # Plot the random fits
    for i in xrange(nsamples):
        # Parameter values at this walker/step
        p = sampler.chain[0, walker_indices[i], step_indices[i], :]
        plot_args = {'color': 'black', 'alpha': 0.1, 'linewidth': 1}
        # Plot a single fit
        gf.plot_func(p, obs_ix=0, ax=ax, plot_args=plot_args,
                     normalize_to_f0=True)

    # Plot the maximum a posteriori fit
    """
    maxp_flat_ix = np.argmax(sampler.lnprobability[0])
    maxp_ix = np.unravel_index(maxp_flat_ix,
                               sampler.lnprobability[0].shape)
    maxp = sampler.lnprobability[0][maxp_ix]
    plot_args = {'color': 'blue', 'alpha': 1, 'linewidth':0.5}
    gf.plot_func(sampler.chain[0, maxp_ix[0], maxp_ix[1]], ax=ax,
                 obs_ix=0, plot_args=plot_args, normalize_to_f0=True)
    """

    plt.subplots_adjust(left=0.14, bottom=0.14, hspace=0.35, wspace=0.5)
    ax.set_title('NBD-%sC-Bax' % residue, fontsize=fontsize)
    if row_ix == (axarr.shape[0] - 1):
        ax.set_xlabel(r'Time (sec $\times 10^3$)')
    if col_ix == 0:
        ax.set_ylabel('NBD F/$F_0$')
    ax.set_xticks(np.linspace(0, 4000, 5))
    ax.set_xticklabels([int(n) for n in np.linspace(0, 4, 5)])
    format_axis(ax)
    ax.set_xlim(-110, 4000)
    # Set ylim depending on the residue
    if residue == '47':
        ax.set_ylim(0.5, 1.)
    elif residue == '120':
        ax.set_yticks(np.linspace(1, 5, 5))
    elif residue == '179':
        ax.set_ylim(1, 2)

plt.savefig('data1_nbd_example_%sconf_fits.pdf' % num_confs)
plt.savefig('data1_nbd_example_%sconf_fits.png' % num_confs)

