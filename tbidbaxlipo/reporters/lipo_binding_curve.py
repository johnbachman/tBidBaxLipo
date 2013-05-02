from tbidbaxlipo.models.one_cpt import Builder
from pylab import *
from pysb.integrate import odesolve

b = Builder()
b.build_model_t()

# Do liposome titration
lipo_concs = logspace (-9, 3, 50)
t = linspace(0, 500, 200)

ion()
figure()
tBid_0 = b.model.parameters['tBid_0'].value
Bax_0 = b.model.parameters['Bax_0'].value

tBid_fracs_bound = np.zeros(len(lipo_concs))
Bax_fracs_bound = np.zeros(len(lipo_concs))

for i, lipo_conc in enumerate(lipo_concs):
    b.model.parameters['Vesicles_0'].value = lipo_conc
    x = odesolve(b.model, t)
    plot(t, x['mtBid'] / tBid_0)
    plot(t, x['mBax'] / Bax_0)
    tBid_fracs_bound[i] = x['mtBid'][-1] / tBid_0
    Bax_fracs_bound[i] = x['mBax'][-1] / Bax_0

show()

figure()
plot(lipo_concs, tBid_fracs_bound)
plot(lipo_concs, Bax_fracs_bound)
show()
