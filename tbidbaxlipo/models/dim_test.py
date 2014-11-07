from tbidbaxlipo.models.one_cpt import Builder
from pysb.integrate import Solver
from pylab import *

bd = Builder()
bd.build_model_nbd_2_conf_dimer()
t = linspace(0, 1e4, 1e3)

sol = Solver(bd.model, t)
sol.run()

ion()
figure()
plot(t, sol.yexpr['NBD'])

figure()
plot(t, (sol.yobs['cBax'] + sol.yobs['mBax_mono'] + sol.yobs['Bax2']) /
        bd.model.parameters['Bax_0'].value)

figure()
bd.model.parameters['Bax_transloc_kf'].value = 1e-4
bd.model.parameters['Bax_mem_dimerization_kf'].value = 1e-2
lipo_concs = np.logspace(np.log10(0.1), np.log10(21.), 20)
for lc in lipo_concs:
    bd.model.parameters['Vesicles_0'].value = lc
    sol.run()
    plot(t, sol.yexpr['NBD'])
