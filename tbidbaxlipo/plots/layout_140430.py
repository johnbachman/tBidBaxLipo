import sys
import os
import tbidbaxlipo.data
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.models.nbd.multiconf import Builder
from tbidbaxlipo.util import fitting

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '140430_cBid85CDAC_BaxNBD_FRET.csv'))

data = np.loadtxt(timecourse_file, delimiter=',')
st = 0
plt.ion()
#plt.figure()
#plt.plot(data[st:,0], data[st:,1])
#plt.plot(data[st:,0], data[st:,2])
#plt.plot(data[st:,4], data[st:,5])
#plt.plot(data[st:,4], data[st:,6])
#plt.plot(data[st:,8], data[st:,9])
#plt.plot(data[st:,8], data[st:,10])

#plt.plot(data[st:,0], data[st:,3])
c68_time = data[:,4]
c68_fret = data[:,7]

#plt.plot(data[st:,8], data[st:,11])
params_dict = {'c0_to_c1_k': 2e-3,
               'c1_scaling': 0.4,
               'c1_to_c2_k': 1e-3,
               'c2_scaling': 0.6,
               'c1_to_c2_k': 1e-3,
               'c3_scaling': 0.5}

builder = Builder(params_dict=params_dict)
builder.build_model_multiconf(5, c68_fret[0], normalized_data=True, reversible=False)

"""
k1 = builder.model.parameters['c0_to_c1_k']
k2 = builder.model.parameters['c1_to_c2_k']
k1_index = builder.model.parameters.index(k1)
k2_index = builder.model.parameters.index(k2)
k1_est_index = builder.estimate_params.index(k1)
k2_est_index = builder.estimate_params.index(k2)
"""

pysb_fit = fitting.fit_pysb_builder(builder, 'NBD', c68_time, c68_fret)
#plt.subplot(1, 2, 2)
plt.figure()
plt.plot(c68_time, c68_fret, linestyle='', marker='.')
plt.plot(c68_time, pysb_fit.ypred)
plt.xlabel('Time (sec)')

