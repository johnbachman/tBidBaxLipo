from tbidbaxlipo.models.one_cpt import Builder
import numpy as np
from matplotlib import pyplot as plt
from pysb.integrate import odesolve

b = Builder()
b.translocate_Bax()
b.basal_Bax_activation()
b.Bax_nmerizes(4)
b.pores_from_Bax_nmers(4)

t = np.linspace(0, 10000, 1000)
x = odesolve(b.model, t)
plt.ion()
plt.figure()
plt.plot(t, x['pores'])
