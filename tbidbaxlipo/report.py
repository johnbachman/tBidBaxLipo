import numpy as np
from pysb.report import reporter

num_samples = 100

@reporter('Inserted Bax')
def inserted_Bax(mcmc_set):
    # Get model from first chain in the set
    model = mcmc_set.chains[0].options.model
    observable_names = [o.name for o in model.observables]
    if 'iBax' not in observable_names:
        return False
    # This will be an empty list if the observable never occurs
    if model.observables['iBax'].species:
        return True
    else:
        return False

@reporter('Bax dimerizes')
def Bax_dimerizes(mcmc_set):
    # Get model from first chain in the set
    model = mcmc_set.chains[0].options.model
    observable_names = [o.name for o in model.observables]
    if 'Bax2' not in observable_names:
        return False
    # This will be an empty list if the observable never occurs
    if model.observables['Bax2'].species:
        return True
    else:
        return False

@reporter('Bax tetramerizes')
def Bax_tetramerizes(mcmc_set):
    # Get model from first chain in the set
    model = mcmc_set.chains[0].options.model
    observable_names = [o.name for o in model.observables]
    if 'Bax4' not in observable_names:
        return False
    # This will be an empty list if the observable never occurs
    if model.observables['Bax4'].species:
        return True
    else:
        return False

@reporter('tBid/Bax increases monotonically')
def tBid_Bax_monotonically_increasing(mcmc_set):
    """ .. todo:: document the basis for this"""
    num_true = 0
    for i in range(num_samples): 
        x  = mcmc_set.get_sample_simulation()
        if monotonic_increasing(x['tBidBax']):
            num_true += 1
    return float(num_true) / float(num_samples)

@reporter('iBax increases monotonically')
def iBax_monotonically_increasing(mcmc_set):
    """ .. todo:: document the basis for this"""
    num_true = 0
    for i in range(num_samples): 
        x  = mcmc_set.get_sample_simulation()
        if monotonic_increasing(x['iBax']):
            num_true += 1
    return float(num_true) / float(num_samples)


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
        

