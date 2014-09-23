from tbidbaxlipo.models.bid_bax_fret import Builder
import numpy as np
from matplotlib import pyplot as plt
from pysb.integrate import Solver

# Constraints considered for parameters
# Initial partitioning to membranes should have a reverse rate (off-rate) of
# 1.6e-3 (Aisha's paper) with an apparent rate of 0.013:
tBid_mem_kr = 0.0016
tBid_mem_kf = 0.013 - 0.0016

# The equilibrium amount of membrane binding should result in, say
# roughly 10% of Bax bound. With the pseudo first order rate above, we can
# calculate what the reverse rate should be explicitly:
def calc_kr(kf, eq): return kf * ((1 - eq)/float(eq))
trans_kf = 0.2
Bax_mem_eq = 0.9
Bax_mem_kr = calc_kr(trans_kf, Bax_mem_eq)

# Bid must undergo a conformational change (but here is preincubated so that
# should not be a factor). Perhaps need to model an initial condition with most
# Bid at membranes but also jumping???

params_dict = {
    'tBid_mBax_kf': 0.001,
    'tBid_mBax_kr': 0.01,
    #'tBidmBax_tBidiBax_kf': 0.001,
    #'tBidmBax_tBidiBax_kr': 0.0001,
    #'tBidiBax_tBid_iBax_kf': 0.01,
    #'tBidiBax_tBid_Bax_kr': 0.001,
    'tBidmBax_tBid_iBax_kf': 0.005,
    'iBax_to_mBax_k': 0.0001,

    'tBid_to_mem_kf': tBid_mem_kf,
    'tBid_to_mem_kr': tBid_mem_kr,
    'Bax_to_mem_kf': trans_kf,
    'Bax_to_mem_kr': Bax_mem_kr,

    'c1_nbd': 20,
    'c2_nbd': 3,
    'c3_nbd': 5,
    'c1_fret': 30,
    'c2_fret': 20,
}

builder = Builder(params_dict=params_dict)

builder.build_model_fret6()

tmax = 4050
t = np.linspace(0, tmax, 10000)
s = Solver(builder.model, t)
s.run()

plt.ion()
plt.close('all')


plt.figure('NBD')
plt.plot(t, s.yexpr['NBD'])
plt.xlim([-50, tmax])
plt.figure('FRET')
plt.plot(t, s.yexpr['FRET'])
plt.xlim([-50, tmax])

plt.figure('Curves')
plt.plot(t, s.yobs['cBax'], label='cBax')
plt.plot(t, s.yobs['mBax_free'], label='mBax_free')
plt.plot(t, s.yobs['iBax_free'], label='iBax_free')
plt.plot(t, s.yobs['mtBid_free'], label='mtBid_free')
plt.plot(t, s.yobs['ctBid'], label='ctBid')
plt.plot(t, s.yobs['tBidmBax'], label='tBidmBax')
plt.plot(t, s.yobs['tBidiBax'], label='tBidiBax')
plt.xlim([-50, tmax])
plt.legend(loc='upper right', prop={'size':8})

def test_Bax_NBD_is_1_without_Bid():
    """In the absence of Bid, or with no binding between them, Bax NBD should
    just be 1."""
    pass

def test_translocation_timescale():
    """Translocation of both tBid and Bax, in the absence of binding between them,
    should be done within a few seconds."""

