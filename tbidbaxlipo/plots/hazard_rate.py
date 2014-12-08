from tbidbaxlipo.plots.layout_140311 import bgsub_norm_wells as data
from tbidbaxlipo.util.plate_assay import TIME, VALUE
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import moving_average, set_fig_params_for_publication, \
        format_axis
from tbidbaxlipo.util import fitting
import sys

plt.ion()
plt.close('all')

well = 'A1'
tc = data['A1']
t = tc[TIME]
y = tc[VALUE]
#plt.plot(t, y)

k1 = fitting.Parameter(1e-4)
def fit_func1(t):
    return (1 - np.exp(-k1()*t))

fmax = fitting.Parameter(0.3)
k2 = fitting.Parameter(1e-4)
def fit_func2(t):
    return fmax() * (1 - np.exp(-k2()*t))

res1 = fitting.fit(fit_func1, [k1], y, t)
res2 = fitting.fit(fit_func2, [fmax, k2], y, t)

# Plot the raw curve and an exponential fit
set_fig_params_for_publication()
plt.figure(figsize=(1.5, 1.5), dpi=300)
plt.plot(t, y, color='k', linewidth=1)
plt.plot(t, fit_func1(t), color='g', linewidth=1)
plt.plot(t, fit_func2(t), color='r', linewidth=1)
plt.xlim([0, 11.3e3])
plt.xlabel(r'Time (sec $\times 10^3$)')
plt.ylabel(r'\% Release')
ax = plt.gca()
format_axis(ax)
ax.set_xticks(np.linspace(0, 1e4, 6))
ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
plt.subplots_adjust(bottom=0.19, left=0.21, top=0.94, right=0.94)
leg = plt.legend(['Data', 'Fit, constant rate', 'Fit, decr. rate'],
                 loc='lower right', prop={'size':6}, frameon=False,
                 handlelength=1.5, handletextpad=0)
#label_padding = 2
#ax.xaxis.labelpad = label_padding
#ax.yaxis.labelpad = label_padding

plt.savefig('fig_fmax_fit_comparison.pdf')

def failure_rate(t, y):
    dt = np.diff(t)
    dy = -np.diff(y)
    fr = dy / (dt * y[1:])
    return fr

plt.figure(figsize=(1.5, 1.5), dpi=300)
#plt.plot(t[1:], failure_rate(t, 1 - y))
plt.plot(t[30:], moving_average(failure_rate(t, 1-y), n=30), color='k',
                linewidth=1)
plt.plot(t[30:], moving_average(failure_rate(t, 1 - fit_func1(t)), n=30),
                color='g', linewidth=1)
plt.plot(t[30:], moving_average(failure_rate(t, 1 - fit_func2(t)), n=30),
                color='r', linewidth=1)
plt.ylabel(r'Hazard rate $(\times 10^{-5})$')
plt.xlabel(r'Time (sec $\times 10^3)$')
ax = plt.gca()
format_axis(ax)
ax.set_xticks(np.linspace(0, 1e4, 6))
ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
ax.set_yticks(np.linspace(-5e-5, 25e-5, 7))
ax.set_yticklabels([int(f) for f in np.linspace(-5, 25, 7)])
#ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])

plt.subplots_adjust(bottom=0.19, left=0.21, top=0.94, right=0.94)
#leg = plt.legend(['Data', 'Fit, constant rate', 'Fit, decr. rate'],
#                 loc='lower right', prop={'size':6}, frameon=False,
#                 handlelength=1.5, handletextpad=0)
#label_padding = 2
#ax.xaxis.labelpad = label_padding
#ax.yaxis.labelpad = label_padding

plt.savefig('fig_hazard_rate.pdf')

