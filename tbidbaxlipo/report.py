import numpy as np
from pysb.report import reporter, Result

num_samples = 100

def write_species(mcmc_set_name, model):
    """Write the species list of the model to a file.

    If the file does not exist, the function does nothing, but returns the name
    of the file.

    Parameters
    ----------
    mcmc_set_name : string
        The name of the MCMC set, which serves as the basename of the file to
        be written out.
    model : pysb.core.Model
        The model whose species are to be written to a file.

    Returns
    -------
    string
        The name of the species list filename.
    """
    species_filename = '%s_species.txt' % mcmc_set_name
    try:
        with open(species_filename, 'w') as f:
            #f.write('\n'.join([str(s) for s in model.rules]))
            #f.write('\n')
            f.write('\n'.join([str(s) for s in model.species]))
    except IOError as e:
        pass
    return species_filename

@reporter('Inserted Bax')
def inserted_Bax(mcmc_set):
    # Get model from first chain in the set
    model = mcmc_set.chains[0].options.model
    observable_names = [o.name for o in model.observables]

    # Write species list to a file
    species_filename = write_species(mcmc_set.name, model)

    if 'iBax' not in observable_names:
        return Result(False, species_filename)
    # This will be an empty list if the observable never occurs
    if model.observables['iBax'].species:
        return Result(True, species_filename)
    else:
        return Result(False, species_filename)

@reporter('Bax dimerizes')
def Bax_dimerizes(mcmc_set):
    # Get model from first chain in the set
    model = mcmc_set.chains[0].options.model
    observable_names = [o.name for o in model.observables]

    # Write species list to a file
    species_filename = write_species(mcmc_set.name, model)

    if 'Bax2' not in observable_names:
        return Result(False, species_filename)
    # This will be an empty list if the observable never occurs
    if model.observables['Bax2'].species:
        return Result(True, species_filename)
    else:
        return Result(False, species_filename)

@reporter('Bax tetramerizes')
def Bax_tetramerizes(mcmc_set):
    # Get model from first chain in the set
    model = mcmc_set.chains[0].options.model
    observable_names = [o.name for o in model.observables]

    # Write species list to a file
    species_filename = write_species(mcmc_set.name, model)

    if 'Bax4' not in observable_names:
        return Result(False, species_filename)
    # This will be an empty list if the observable never occurs
    if model.observables['Bax4'].species:
        return Result(True, species_filename)
    else:
        return Result(False, species_filename)

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

