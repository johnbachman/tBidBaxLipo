from pysb import *
import nbd_model_shared
import numpy as np

Model()

nbd_model_shared.declare_shared_components()

Parameter('c3_insertion_rate', 0.004053)
Parameter('c62_insertion_rate', 0.00112)
Parameter('c120_insertion_rate', 0.001693)
Parameter('c122_insertion_rate', 0.001227)
Parameter('c126_insertion_rate', 0.002)

Rule('c3_insertion', Bax(c3='s') >> Bax(c3='m'), c3_insertion_rate)
Rule('c62_insertion', Bax(c62='s') >> Bax(c62='m'), c62_insertion_rate)
Rule('c120_insertion', Bax(c120='s') >> Bax(c120='m'), c120_insertion_rate)
Rule('c122_insertion', Bax(c122='s') >> Bax(c122='m'), c122_insertion_rate)
Rule('c126_insertion', Bax(c126='s') >> Bax(c126='m'), c126_insertion_rate)

def random_initial_values(num_sets):
    """Generate a random sample of initial values for parameter estimation.

    Calls nbd_model_shared.random_initial_values to get randomized initial
    values for the scaling factor parameters, and then concatenates
    randomized initial values for the transition rate parameters. 5 scaling
    factor parameters + 5 transition rate parameters = 10 parameter total.

    Logs of the transition rate parameters are sampled from a uniform
    distribution ranging from [-6, 0]. These values are exponentiated to
    return values ranging from 10^-6 to 1.

    The argument num_sets specifies the number of sets of parameter values
    to produce. For example, to run 10 different MCMC chains, generate 10 
    parameter sets.
    """

    lower_bound = -6
    upper_bound = 0
    initial_values_list = []
    for i in range(0, num_sets):
        initial_values_list.append(10 ** np.random.uniform(low=lower_bound,
                                                     high=upper_bound, size=5))

    # Concatenate with the randomized initial values from nbd_model_shared
    return np.concatenate((nbd_model_shared.random_initial_values(num_sets),
                           initial_values_list), axis=1)

