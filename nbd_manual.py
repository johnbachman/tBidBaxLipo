from tBid_Bax_1c import tBid_Bax_1c

import nbd_analysis as nbd
import matplotlib.pyplot as plt
from pysb import *
from pysb.integrate import odesolve

params_dict = {
'tBid_transloc_kf':1e-1,
'tBid_transloc_kr':0,
'Bax_transloc_kf':1e-1,
'Bax_transloc_kr':100,
'tBid_mBax_kf':1,
'tBid_mBax_kr':0.01,
'tBid_iBax_kc':10,
'iBax_reverse_k':1e-1,
'Bax_dimerization_kf':1e-3,
'Bax_dimerization_kr':1e-2,
'Bax_tetramerization_kf':1e-3,
'Bax_tetramerization_kr':1e-4
}

# After 2000 steps of fitting: very good fit!
"""
[(Parameter(name='tBid_transloc_kf', value=0.5), 0.046433226280901206),
 (Parameter(name='tBid_transloc_kr', value=0), 0.0),
 (Parameter(name='Bax_transloc_kf', value=0.5), 5.007883980982208),
 (Parameter(name='Bax_transloc_kr', value=100), 171.1669217159442),

 (Parameter(name='tBid_mBax_kf', value=0.20000000000000001),
  0.062823108082142434),
 (Parameter(name='tBid_mBax_kr', value=0.01), 0.036115738840543629),
 (Parameter(name='tBid_iBax_kc', value=10), 11.621431266275815),
 (Parameter(name='iBax_reverse_k', value=0.10000000000000001),
  0.2713605496448912),

 (Parameter(name='Bax_dimerization_kf', value=0.001), 0.00039997477906346336),
 (Parameter(name='Bax_dimerization_kr', value=0.01), 0.013749131382056693),
 (Parameter(name='Bax_tetramerization_kf', value=0.001),
  0.00093066359281583318),
 (Parameter(name='Bax_tetramerization_kr', value=0.0001),
  0.0002291780988903229)]
"""

# by name:
[('tBid_transloc_kf', 0.046433226280901206),
 ('tBid_transloc_kr', 0.0),
 ('Bax_transloc_kf', 5.007883980982208),
 ('Bax_transloc_kr', 171.1669217159442),
 ('tBid_mBax_kf', 0.062823108082142434),
 ('tBid_mBax_kr', 0.036115738840543629),
 ('tBid_iBax_kc', 11.621431266275815),
 ('iBax_reverse_k', 0.2713605496448912),
 ('Bax_dimerization_kf', 0.00039997477906346336),
 ('Bax_dimerization_kr', 0.013749131382056693),
 ('Bax_tetramerization_kf', 0.00093066359281583318),
 ('Bax_tetramerization_kr', 0.0002291780988903229)]



# f/r 1e6 M^-1 = 1 / 1e6 M = 1uM kD

m1c = tBid_Bax_1c(params_dict=params_dict)
m1c.build_model2()
model = m1c.model

tspan = nbd.time_c62
nbd_avgs, nbd_stds = nbd.calc_norm_avg_std()
ydata_norm = nbd_avgs[1] # NBD_62c 

# Plot "Before" curves
plt.plot(tspan, ydata_norm)
before = odesolve(m1c.model, tspan)
#before_array = before.view().reshape(len(tspan), len(before.dtype))
#iBax_before = 2*before_array[:,6]
#iBax_before_norm = iBax_before / max(iBax_before)
Baxc62 = (2*before['Bax2'] + 4*before['Bax4']) / 100
Baxc62_norm = Baxc62 / max(Baxc62)
#Baxc62 = 
plt.plot(tspan, before['iBax']/100, label='iBax')
plt.plot(tspan, before['mBax']/100, label='mBax')
plt.plot(tspan, before['Bax2']/100, label='Bax2')
plt.plot(tspan, before['Bax4']/100, label='Bax4')
plt.plot(tspan, Baxc62, label='c62')
plt.plot(tspan, Baxc62_norm, label='c62_norm')
plt.legend(loc='lower right')

"""
The idea that I'm going for is that if iBax accumulates quickly,
this will be followed by rapid increase in Bax dimers. The rate of Bax dimer
increase is  fastest when iBax is highest, but then it starts to slow down
as iBax starts to go down (iBax is being consumed by Bax2).

However, if we set iBax 
If we set the parameters such that dimerization happens quickly,
"""
