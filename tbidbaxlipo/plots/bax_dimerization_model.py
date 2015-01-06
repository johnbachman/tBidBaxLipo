from tbidbaxlipo.models.one_cpt import Builder
from pysb.integrate import Solver
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.util import fitting

b = Builder()
#b.build_model_nbd_2_conf()
b.build_model_nbd_2_conf_dimer()

t = np.linspace(0, 100000, 1e3)
sol = Solver(b.model, t)

plt.ion()
plt.close('all')

bax_concs = np.logspace(0, 3, 30)
k_list = []
fmax_list = []

for bax_ix, bax_conc in enumerate(bax_concs):
    b.model.parameters['Bax_0'].value = bax_conc
    #b.model.parameters['Vesicles_0'].value = bax_conc
    sol.run()
    plt.figure('NBD')
    plt.plot(t, sol.yexpr['NBD'])
    fmax = fitting.Parameter(5.)
    k = fitting.Parameter(5e-4)
    def exp_func(t):
        return fmax()*(1 - np.exp(-k()*t)) + 1.
    fitting.fit(exp_func, [k, fmax], sol.yexpr['NBD'], t)
    k_list.append(k())
    fmax_list.append(fmax())

plt.figure('k')
plt.plot(bax_concs, k_list)
ax = plt.gca()
ax.set_xscale('log')

plt.figure('fmax')
plt.plot(bax_concs, fmax_list)
ax = plt.gca()
ax.set_xscale('log')


