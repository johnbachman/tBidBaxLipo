import nbd_analysis as nbd
import matplotlib.pyplot as plt
from pysb import *
from pysb.integrate import odesolve
from nbd_model import model

params_dict = { }

#tspan = nbd.time_c62
nbd_avgs, nbd_stds = nbd.calc_norm_avg_std()
#ydata_norm = nbd_avgs[1] # NBD_62c 

# Plot data curves
plt.ion()
plt.plot(nbd.time_other, nbd_avgs[0], 'r,', label='c3 data', alpha=0.2)
plt.plot(nbd.time_c62, nbd_avgs[1], 'g,', label='c62 data', alpha=0.2)
plt.plot(nbd.time_other, nbd_avgs[2], 'b,', label='c120 data', alpha=0.2)
plt.plot(nbd.time_other, nbd_avgs[3], 'm,', label='c122 data', alpha=0.2)
plt.plot(nbd.time_other, nbd_avgs[4], 'k,', label='c126 data', alpha=0.2)

# Plot model curves
tspan = nbd.time_other
x = odesolve(model, tspan)
plt.plot(tspan, x['Baxc3'] * model.parameters['c3_scaling'].value,
         'r', label='c3 model')
plt.plot(tspan, x['Baxc62'], 'g', label='c62 model')
plt.plot(tspan, x['Baxc120'], 'b', label='c120 model')
plt.plot(tspan, x['Baxc122'], 'm', label='c122 model')
plt.plot(tspan, x['Baxc126'], 'k', label='c126 model')
plt.show()
plt.legend(loc='lower right')

