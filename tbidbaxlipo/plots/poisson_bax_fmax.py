from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import poisson
from tbidbaxlipo.util import set_fig_params_for_publication
from tbidbaxlipo.plots.layout_140311 import get_twoexp_fmax_arr
from tbidbaxlipo.util import fitting
from matplotlib import ticker

plt.ion()
set_fig_params_for_publication()

# First, plot the figure showing the relationship between the Poisson
# distribution and the percent permeabilized
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
ax.set_yticks([0, 0.5, 1.0])

plt.ylim([-0.08, 1.05])
plt.xlabel('Bax*/Lipo')
plt.ylabel('Max release')

plt.savefig('poisson_bax_fmax.pdf')

# Now, plot best fit of 140311 Fmax curve with Poisson funcs
(fmax_arr, conc_list) = get_twoexp_fmax_arr()
fmax_means = np.mean(fmax_arr, axis=0)
bax_ratios = conc_list / 5.16
log_ratios = np.log10(bax_ratios)
fmax_arr[:, 0] = [0, 0, 0]

plt.figure(figsize=(1.5, 1.5), dpi=300)
plt.subplots_adjust(left=0.21, bottom=0.17, top=0.89, right=0.96)
line_data = plt.errorbar(bax_ratios, fmax_means, yerr=np.std(fmax_arr, axis=0))
plt.xlabel('[Bax]/[Liposome]')
plt.ylabel('Predicted max release')
plt.xlim([1, 150])
plt.ylim([-0.015, 1.015])
ax = plt.gca()
ax.set_xscale('log')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.yaxis.set_tick_params(direction='out')
ax.xaxis.set_tick_params(which='both', direction='out')
ax.xaxis.labelpad = label_padding
ax.yaxis.labelpad = label_padding
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
# Because the poisson.sf only differs at discrete values of the pore size,
# leastsq fails. So instead we'll just iterate over a range of pore sizes and
# find the one with the minimum error.
def poisson_pore_fit(min_pore_size):
    return np.sum((fmax_means[1:] -
                   poisson.sf(min_pore_size, bax_ratios[1:])) ** 2)

best_pore_size = 1
best_fit_err = np.inf
for pore_size in range(1, 200):
    err = poisson_pore_fit(pore_size)
    if err < best_fit_err:
        best_fit_err = err
        best_pore_size = pore_size
smooth_bax_ratios= np.logspace(log_ratios[1], log_ratios[-1], 50)
line_4 = plt.plot(smooth_bax_ratios, poisson.sf(4, smooth_bax_ratios),
         color='g')
line_best = plt.plot(smooth_bax_ratios,
                     poisson.sf(best_pore_size, smooth_bax_ratios),
                     color='r')

leg = ax.legend(['Data', 'Min 4', 'Min 33 (best fit)'], loc=1, ncol=3,
          bbox_to_anchor=(0.0, 0.9, 1, 0.1),
          bbox_transform=plt.gcf().transFigure,
          mode='expand', borderpad=0,
          handlelength=1.5, handletextpad=0, columnspacing=0.5,
          frameon=False, prop={'size': 6})

plt.savefig('poisson_bax_fmax_fit.pdf')

"""
ax = plt.axes([0.05, 0.9, 0.9, 0.1])
#plt.subplots_adjust(left=0, top=0, bottom=0, right=0)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
#plt.setp(ax, xticks=[], yticks=[])
plt.hlines(0.5, 0, 0, color='b')
plt.text(0, 0, 'Data')
"""
