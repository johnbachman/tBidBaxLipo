from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import poisson
from tbidbaxlipo.util import set_fig_params_for_publication

plt.ion()
set_fig_params_for_publication()
plt.figure(figsize=(1.5, 1.5), dpi=300)

label_padding = 2
min_pore_size = 4
bax_ratios = [1, 4, 6, 10]
index = np.arange(21)
maxval = -1
for br in bax_ratios:
    plt.plot(index, poisson.pmf(index, br))
    poisson_max = np.max(poisson.pmf(index, br))
    if poisson_max > maxval:
        maxval = poisson_max
plt.xlabel('Bax*/Liposome')
plt.ylabel('Probability')
ax = plt.gca()
line_top = ax.get_ylim()[1]
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('none')
ax.set_yticklabels([])
ax.xaxis.labelpad = label_padding
ax.yaxis.labelpad = 0
ax.yaxis.set_tick_params(direction='out')
ax.xaxis.set_tick_params(direction='out')
plt.vlines(4, 0, line_top, linestyle='--')
plt.ylim([0, line_top])
plt.subplots_adjust(left=0.12, bottom=0.17)

ax = plt.axes([0.55, 0.55, 0.3, 0.3])
bax_ratios2 = np.linspace(0, 20, 50)
plt.plot(bax_ratios2, poisson.sf(4, bax_ratios2), color='k')
for br in bax_ratios:
    plt.plot(br, poisson.sf(4, br), marker='o', markersize=4)
ax.xaxis.labelpad = label_padding
ax.yaxis.labelpad = label_padding
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.yaxis.set_tick_params(direction='out')
ax.xaxis.set_tick_params(direction='out')
plt.ylim([-0.08, 1.05])
plt.xlabel('Bax*/Lipo')
plt.ylabel('Max release')
