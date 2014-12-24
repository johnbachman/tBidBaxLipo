from tbidbaxlipo.models.one_cpt import Builder
from pysb.integrate import Solver
import numpy as np
from matplotlib import pyplot as plt

b = Builder()
b.build_model_nbd_2_conf()
#b.build_model_nbd_2_conf_dimer()

t = np.linspace(0, 10000, 1e3)
sol = Solver(b.model, t)

plt.ion()
plt.close('all')

bax_concs = np.logspace(0, 3, 30)
for bax_conc in bax_concs:
    b.model.parameters['Bax_0'].value = bax_conc
    sol.run()
    plt.figure('NBD')
    plt.plot(t, sol.yexpr['NBD'])
