"""
A script to parse the NBD fluorescence data from the Excel spreadsheet sent by
Justin Kale into Python data files.
"""

import numpy as np
from tbidbaxlipo.util import fitting
from matplotlib import pyplot as plt
from pysb import *
from pysb.integrate import Solver

# Define a simple model
Model()
Monomer('Bid', ['loc'], {'loc': ['c', 'm', 'i']})
Initial(Bid(loc='c'), Parameter('Bid_0', 1))
Rule('Bid_binds_liposomes', Bid(loc='c') >> Bid(loc='m'),
        Parameter('Bid_c_to_m_k', 0.215))
Rule('Bid_conf_change', Bid(loc='m') >> Bid(loc='i'),
        Parameter('Bid_m_to_i_k', 0.215))
Observable('Bid_m', Bid(loc='m'))
Observable('Bid_i', Bid(loc='i'))

# Fit the data
data = np.loadtxt('cBid1uM.txt')
time = data[:2500,0]
y = data[:2500,1]

k1 = fitting.Parameter(0.001)
fmax1 = fitting.Parameter(0.3)
k2 = fitting.Parameter(0.001)
fmax2 = fitting.Parameter(0.3)

def single_exp(t):  return ((fmax1()*(1 - np.exp(-k1()*t))))
def double_exp(t):  return ((fmax1()*(1 - np.exp(-k1()*t)))  +
                            (fmax2()*(1 - np.exp(-k2()*t))))

solver = Solver(model, time)

def pysb_model(t):
    #scaling_Bid_m = fmax1()
    Bid_c_to_m_k.value = k1()
    #scaling_Bid_i = fmax2()
    Bid_m_to_i_k.value = k2()
    solver.run()
    return (fmax1() * solver.yobs['Bid_m']) + (fmax2() * solver.yobs['Bid_i'])

fitting.fit(pysb_model, [fmax1, k1, fmax2, k2], y, time)

print "k1 = %f" % k1()
print "fmax1 = %f" % fmax1()
print "k2 = %f" % k2()
print "fmax2 = %f" % fmax2()

plt.ion()
plt.figure()
plt.plot(time, y)
plt.plot(time, pysb_model(time), 'r')
plt.show()


