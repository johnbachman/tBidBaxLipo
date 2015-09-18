from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis, \
                             fontsize

set_fig_params_for_publication()

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
        ('179', 1),
        ]
# Get data for NBD-54C-Bax

end_ix = -1
fig, axarr = plt.subplots(2, 3, sharex=True, figsize=(3, 2), dpi=300)
for plot_index, plot_tuple in enumerate(curves_to_plot):
    residue, rep = plot_tuple
    (row_ix, col_ix) = np.unravel_index(plot_index, axarr.shape)
    ax = axarr[row_ix, col_ix]
    t = df[('Bid', 'NBD', residue, rep, 'TIME')][:end_ix]
    v = df[('Bid', 'NBD', residue, rep, 'VALUE')][:end_ix]
    norm_v = v / float(v[0])
    ax.plot(t, norm_v, label=residue)
    ax.set_xlim(0, 4000)
    plt.subplots_adjust(left=0.14, bottom=0.14, hspace=0.35, wspace=0.5)
    ax.set_title('NBD-%sC-Bax' % residue, fontsize=fontsize)
    if row_ix == (axarr.shape[0] - 1):
        ax.set_xlabel('Time (sec)')
    if col_ix == 0:
        ax.set_ylabel('NBD F/$F_0$')
    ax.set_xticks(np.linspace(0, 4000, 5))
    ax.set_xticklabels([int(n) for n in np.linspace(0, 4, 5)])
    format_axis(ax)

#ax.set_xlabel('Time (sec)')
#ax.set_ylabel('NBD $F/F_0$')
#ax.set_xticks([0, 1000, 2000, 3000, 4000])
#ax.set_yticks([1, 2, 3, 4, 5])
#format_axis(ax)
#leg = plt.legend(loc='right', fontsize=6)
#leg.draw_frame(False)
plt.savefig('data1_nbd_example_curves.pdf')
plt.savefig('data1_nbd_example_curves.png')

