import numpy as np
from pysb.report import reporter

@reporter('tBid/Bax increases monotonically')
def tBid_Bax_monotonically_increasing(mcmc):
    """ .. todo:: document the basis for this"""
    x = mcmc.simulate(position=mcmc.positions[-1], observables=True)
    return monotonic_increasing(x['tBidBax'])

@reporter('iBax increases monotonically')
def iBax_monotonically_increasing(mcmc):
    """ .. todo:: document the basis for this"""
    x = mcmc.simulate(position=mcmc.positions[-1], observables=True)
    return monotonic_increasing(x['iBax'])

# Helper functions
# ================
def monotonic_increasing(x):
    # TODO rewrite so doesn't allow fixed, unchanging values
    dx = np.diff(x)
    return bool(np.all(dx >= 0))
           
def monotonic_decreasing(x):
    # TODO rewrite so doesn't allow fixed, unchanging values
    dx = np.diff(x)
    return bool(np.all(dx <= 0))
        
