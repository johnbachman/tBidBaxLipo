from tbidbaxlipo.models.one_cpt import Builder
from pysb.integrate import Solver
import numpy as np
from matplotlib import pyplot as plt
from preprocess_data import time_126, data_126

"""
Parameter('Vesicles_0', 5.0),
 Parameter('tBid_0', 20.0),
 Parameter('Bax_0', 100.0),
 Parameter('Bax_transloc_kf', 0.01),
 Parameter('Bax_transloc_kr', 0.1),
 Parameter('basal_Bax_kf', 0.002),
 Parameter('Bax_ins_dimerization_kf', 0.002),
 Parameter('c0_scaling', 1.0),
 Parameter('c1_scaling', 4.0),
 Parameter('bax_fret_scaling', 50.0),
 ])
"""

params_dict = {
        'c1_scaling': 4.8,
        'bax_fret_scaling': 60.,
        'basal_Bax_kf': 0.004,
        'Bax_DAC_0': 20.0,
        'Bax_NBD_0': 100.0,
        'Bax_0': 0.0,
        }

bd = Builder(params_dict=params_dict)

d = {'baxtranslocation': 1,
     'activation': 1,
     'dimerization': 1,
     'nbd': 2,
     'baxfret': 1,
     }
bd.build_model_from_dict(d)

#t = np.linspace(0, 1e4, 1e3)
t = time_126
sol = Solver(bd.model, t)
sol.run()

plt.ion()
plt.close('all')
plt.figure(figsize=(4, 8))

plt.subplot(2, 1, 1)
plt.plot(t, sol.yexpr['NBD'])
plt.plot(time_126, data_126[0, 0, :])

plt.subplot(2, 1, 2)
plt.plot(t, sol.yexpr['BaxFRET'])
plt.plot(time_126, data_126[0, 1, :])

