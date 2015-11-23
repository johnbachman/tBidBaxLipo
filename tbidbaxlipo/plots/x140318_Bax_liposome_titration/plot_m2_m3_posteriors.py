import cPickle
import sys
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis

#filelist = sys.argv[1:]

set_fig_params_for_publication()
plt.ion()

filename = 'mcmc/pt_140318_1c_Baxtr1Activ1Rever1Nbd1.mcmc'
print "Loading %s" % filename
with open(filename) as f:
    (gf, sampler) = cPickle.load(f)
ibax_reverse_k = gf.builder.model.parameters['iBax_reverse_k']
ibax_index = gf.builder.estimate_params.index(ibax_reverse_k)
# Get samples
ibax_samples = sampler.flatchain[0, :, ibax_index]
ibax_mean = np.mean(ibax_samples)
ibax_sd = np.std(ibax_samples, ddof=1)

plt.figure(figsize=(1.1, 1), dpi=300)
plt.hist(ibax_samples, bins=20, normed=True)
plt.xlabel('log$_{10}(k_{r2})$')
plt.ylabel('Density')
plt.xticks([-12, -10, -8, -6, -4])
plt.yticks([0, 0.1, 0.2, 0.3, 0.4])
plt.subplots_adjust(bottom=0.25, left=0.28)
ax = plt.gca()
format_axis(ax)
plt.savefig('%s.ibax_reverse.pdf' % filename)

filename = 'mcmc/pt_140318_ls_Baxtr1Activ1Nbd1.mcmc'
print "Loading %s" % filename
with open(filename) as f:
    (gf, sampler) = cPickle.load(f)
nsites_k = gf.builder.model.parameters['sites_per_liposome']
nsites_index = gf.builder.estimate_params.index(nsites_k)
# Get samples
nsites_samples = sampler.flatchain[0, :, nsites_index]
nsites_mean = np.mean(nsites_samples)
nsites_sd = np.std(nsites_samples, ddof=1)

plt.figure(figsize=(1.1, 1), dpi=300)
plt.hist(nsites_samples, range=(1, 3), bins=20, normed=True)
plt.xlabel('log$_{10}$(sites/lipo)')
plt.ylabel('Density')
#plt.xticks([-12, -10, -8, -6, -4])
#plt.yticks([0, 0.1, 0.2, 0.3, 0.4])
plt.subplots_adjust(bottom=0.25, left=0.28)
ax = plt.gca()
format_axis(ax)
plt.savefig('%s.nsites.pdf' % filename)

