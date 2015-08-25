from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis

#plt.ion()

fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
ax = fig.gca()


curves_to_plot = [
        ('3', 1),
        #('5', 1),
        #('15', 1),
        ('54', 1),
        #('62', 2),
        ('68', 1),
        ('120', 1),
        ('179', 1),
        #('47', 1),
        ]
# Get data for NBD-54C-Bax

end_ix = -1
for residue, rep in curves_to_plot:
    t = df[('Bid', 'NBD', residue, rep, 'TIME')][:end_ix]
    v = df[('Bid', 'NBD', residue, rep, 'VALUE')][:end_ix]
    norm_v = v / float(v[0])
    ax.plot(t, norm_v, label=residue)

plt.subplots_adjust(left=0.19, bottom=0.19)
ax.set_xlabel('Time (sec)')
ax.set_ylabel('NBD $F/F_0$')
ax.set_xticks([0, 1000, 2000, 3000, 4000])
ax.set_yticks([1, 2, 3, 4, 5])
format_axis(ax)
leg = plt.legend(loc='right', fontsize=6)
leg.draw_frame(False)
plt.savefig('data1_nbd_example_curves.pdf')
plt.savefig('data1_nbd_example_curves.png')
