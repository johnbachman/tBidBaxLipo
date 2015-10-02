from tbidbaxlipo.util import set_fig_params_for_publication, format_axis, \
                             fontsize
from tbidbaxlipo.plots.bid_bim_fret_nbd_release.preprocess_data import \
        df_pre as df
import numpy as np
from matplotlib import pyplot as plt

set_fig_params_for_publication()

curves_to_plot = [('Bid', '54'),
                  ('Bid', '126'),
                  ('Bim', '54'),
                  ('Bim', '126'),]

plt.ion()

fig, axarr = plt.subplots(2, 2, sharex=True, figsize=(3, 3), dpi=300)

# Four subplots
for plot_ix, curve_info in enumerate(curves_to_plot):
    (row_ix, col_ix) = np.unravel_index(plot_ix, axarr.shape)
    act = curve_info[0]
    res = curve_info[1]
    ax1 = axarr[row_ix, col_ix]
    ax2 = ax1.twinx()
    # Get data
    nbd_t = df[(act, 'NBD', res, 1, 'TIME')].values
    nbd_v = df[(act, 'NBD', res, 1, 'VALUE')].values
    fret_t = df[(act, 'FRET', res, 1, 'TIME')].values
    fret_v = df[(act, 'FRET', res, 1, 'VALUE')].values
    rel_t = df[(act, 'Release', res, 1, 'TIME')].values
    rel_v = df[(act, 'Release', res, 1, 'VALUE')].values
    # Set title
    act_name = 'cBid' if act == 'Bid' else act
    ax1.set_title('NBD-%sC-Bax, DAC-%s' % (res, act_name), fontsize=fontsize)
    if row_ix == (axarr.shape[0] - 1):
        ax1.set_xlabel(r'Time (sec $\times 10^{-3}$)')
    if col_ix == 0:
        ax1.set_ylabel('NBD F/$F_0$')
        ax2.set_yticks([1, 1.1, 1.2, 1.3])
        ax2.set_ylim(1, 1.32)
    else:
        ax2.set_yticks([1, 2, 3, 4, 5])
        ax2.set_ylim(1, 5.1)
    ax1.set_xticks(np.linspace(0, 4000, 5))
    ax1.set_xticklabels([int(n) for n in np.linspace(0, 4, 5)])
    ax1.set_xlim(-50, 4000)
    # Label axes
    ax1.set_ylabel('\% FRET, \% Release')
    ax2.set_ylabel('NBD $F/F_0$')
    ax1.set_ylim(0, 70)
    # Plot
    ax1.plot(fret_t, fret_v, color='b')
    ax1.plot(rel_t, rel_v, color='r')
    ax2.plot(nbd_t, nbd_v, color='g')
    # Format axes
    format_axis(ax1)
    format_axis(ax2, yticks_position='right')
    plt.subplots_adjust(left=0.14, bottom=0.14, hspace=0.35, wspace=0.7)

fig.savefig('data2_example_curves.pdf')
fig.savefig('data2_example_curves.png', dpi=300)
