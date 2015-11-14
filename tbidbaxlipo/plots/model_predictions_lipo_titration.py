from tbidbaxlipo.models import one_cpt, lipo_sites
from matplotlib import pyplot as plt
import numpy as np
from pysb.integrate import Solver
from tbidbaxlipo.util import fontsize, set_fig_params_for_publication

set_fig_params_for_publication()

plt.ion()


color_list = [
    '#d73027',
    '#fc8d59',
    '#fee090',
    '#ffffbf',
    '#e0f3f8',
    '#91bfdb',
    '#4575b4'
]

# A function for simulating a titration
def lipo_titration(builder, ax):
    t = np.linspace(0, 1e4, 1e3)
    sol = Solver(bd.model, t)
    lipo_concs = np.logspace(1.5, -1, 7)
    for i, lipo_conc in enumerate(lipo_concs):
        bd['Vesicles_0'].value = lipo_conc
        sol.run()
        ax.plot(t, sol.yexpr['NBD'], color=color_list[i])



# row and column sharing
f, ((ax1, ax2)) = plt.subplots(1, 2, sharex=True, sharey=True,
                                           figsize=(2, 1.1), dpi=300)
ax1.set_ylim(0, 5)
ax2.set_xlim(-200, 1e4)

# Simplest model (M1)
model_dict = {
    'baxtranslocation': 1,
    'activation': 1,
    'nbd': 1
}
params_dict = {
    'Bax_transloc_kf': 1e-3,
    'Bax_transloc_kr': 1e-3,
    'basal_Bax_kf': 1e-3,
    'Vesicles_0': 2,
    'Bax_NBD_0': 100,
    'Bax_0': 0,
    'c1_scaling': 4.8,
}
bd = one_cpt.Builder(params_dict=params_dict)
bd.build_model_from_dict(model_dict)
lipo_titration(bd, ax1)

# Lipo sites model (M3)
model_dict = {
    'builder': 'lipo_sites',
    'baxtranslocation': 1,
    'activation': 1,
    'nbd': 1
}
params_dict = {
    'Bax_transloc_kf': 1e-3,
    'Bax_transloc_kr': 1e-3,
    'basal_Bax_kf': 1e-3,
    'Vesicles_0': 2,
    'Bax_NBD_0': 100,
    'Bax_0': 0,
    'c1_scaling': 4.8,
    'sites_per_liposome': 100,
}
bd = lipo_sites.Builder(params_dict=params_dict)
bd.build_model_from_dict(model_dict)
lipo_titration(bd, ax2)

# Remove ticks
for ax in [ax1, ax2]:
    ax.set_xticks([])
    ax.set_yticks([])

f.subplots_adjust(hspace=0, wspace=0, bottom=0.18)

# Set common labels
f.text(0.5, 0.06, 'Time', ha='center', va='center', fontsize=fontsize)
f.text(0.06, 0.5, r'\% Inserted Bax', ha='center', va='center',
         rotation='vertical', fontsize=fontsize)

vert = 0.5
horiz = 5000
ax1.text(horiz, vert, 'Partitioning (M1)', ha='center',
         fontsize=fontsize)
ax2.text(horiz, vert, 'Lipo sites (M3)', fontsize=fontsize,
         ha='center')

f.savefig('model_predictions_lipo_titration.pdf')
f.savefig('model_predictions_lipo_titration.png')

