from pysb import *
from pysb.util import alias_model_components
import numpy as np

# The baseline state--all sites are soluble
sites_initial_state = {
        'c3':'s',
       'c62':'s',
      'c120':'s',
      'c122':'s',
      'c126':'s',
      'c184':'s'}

def declare_shared_components():

    Monomer('Bax', ['c3', 'c62', 'c120', 'c122', 'c126', 'c184'],
        {'c3': ['s', 'm'],
         'c62': ['s', 'm'],
         'c120': ['s', 'm'],
         'c122': ['s', 'm'],
         'c126': ['s', 'm'],
         'c184': ['s', 'm']})

    Parameter('Bax_0', 1)
    Parameter('c3_scaling', 0.7584)
    Parameter('c62_scaling', 0.9204)
    Parameter('c120_scaling', 0.975)
    Parameter('c122_scaling', 0.952)
    Parameter('c126_scaling', 0.966)

    alias_model_components()

    Observable('Baxc3', Bax(c3='m'))
    Observable('Baxc62', Bax(c62='m'))
    Observable('Baxc120', Bax(c120='m'))
    Observable('Baxc122', Bax(c122='m'))
    Observable('Baxc126', Bax(c126='m'))
    Observable('Baxc184', Bax(c184='m'))

    Initial(Bax(**sites_initial_state), Bax_0)

num_parameters = 5

def random_initial_values(num_sets):
    """Get a matrix of random initial values for the scaling parameters.
    
    The argument num_sets specifies the number of sets of parameter sets
    to produce. For example, to run 10 different MCMC chains, generate 10
    parameter sets.

    Note that because the initial value parameter Bax_0 is not estimated,
    it does not have an included random initial value.
    """
    initial_values_list = []

    for i in range(0, num_sets):
        initial_values_list.append(np.random.uniform(low=0.5, high=1.1, size=num_parameters))

    return initial_values_list


