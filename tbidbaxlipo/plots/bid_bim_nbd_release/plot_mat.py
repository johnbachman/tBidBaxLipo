import pickle
from matplotlib import pyplot as plt
from tbidbaxlipo.util import format_axis, set_fig_params_for_publication
import numpy as np

with open('density_mx.pck') as f:
    density_mx = pickle.load(f)

ncols = 19

plt.ion()

residues = ['3', '5', '15', '54', '62']
fig = plt.figure(figsize=(2, 2), dpi=300)
ax = fig.gca()
ax.matshow(density_mx.T, cmap='gray', extent=(-4, -1, 0, 19), aspect='auto')
format_axis(ax)
lines = np.arange(3, ncols, 4) + 0.5
plt.hlines(lines, -4, -1)
ytick_positions = np.arange(1, 19, 4) + 0.5
ax.set_yticks(ytick_positions)
assert len(residues) == len(ytick_positions)
ax.set_yticklabels(residues)
xtick_positions = np.arange(-4, -0.9, 0.5)
ax.set_xticks(xtick_positions)
plt.show()

