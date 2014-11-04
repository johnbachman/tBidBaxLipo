import pickle
from pylab import *

(gf, samp) = pickle.load(open('pt_140318_nbd_2_conf_3.pck'))

c = samp.chain
ntemps = samp.chain.shape[0]

lkl_m = []
lkl_sd = []

for t_ix in range(ntemps):
    lkl_m.append(np.mean(samp.lnlikelihood[t_ix]))
    lkl_sd.append(np.std(samp.lnlikelihood[t_ix]))

lkl_m = array(lkl_m)
lkl_sd = array(lkl_sd)

ion()
figure()
errorbar(np.log10(samp.betas), lkl_m, yerr=lkl_sd)

figure()
plot(np.log10(samp.betas), np.log10(-lkl_m), marker='o')

#figure()
#plot(samp.lnlikelihood[15,:,:].T, alpha=0.2)
