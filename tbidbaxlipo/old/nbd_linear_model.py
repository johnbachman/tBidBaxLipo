from pysb import *
from pysb.macros import equilibrate
import nbd_model_shared
import numpy as np
from itertools import chain

site_list = ['c3', 'c62', 'c126', 'c122', 'c120']

Model()

nbd_model_shared.declare_shared_components()

def linear_pathway_from_ordering(site_ordering, parameters):
    """Builds a model in which each residue transitions in a given order.

    Sites go from 's' to 'm' in a stepwise fashion, that is, a site cannot
    make the transition unless all previous sites in the order have already
    made the transition. A series of reversible rules are generated
    specifying this sequence based on the given ordering.

    Parameters
    ----------
    site_ordering : list of strings (sites)
        A list of site names (e.g., ['c62', 'c3', ...]) that specify the order
        in which the sites transition from solvated to desolvated (e.g.,
        c62 is first, then c3, etc.).
    """

    # "lhs" denotes the "left hand side" of the rule, that is, the reactant
    # ("rhs" will denote the right hand side, i.e., the product)
    site_states_lhs = nbd_model_shared.sites_initial_state.copy()

    for i, current_site in enumerate(site_ordering):
        site_states_rhs = site_states_lhs.copy()
        site_states_rhs[current_site] = 'm'

        # The equilibration rule defining the insertion step
        equilibrate(Bax(**site_states_lhs), Bax(**site_states_rhs),
                parameters[(2*i):(2*i+2)])   

        # The product of this step becomes the reactant for the next step:
        site_states_lhs = site_states_rhs.copy()

num_residues = 5

## Build the default model # FIXME this should be parameterized
linear_pathway_from_ordering(site_list, [1e-2, 1e-4] * num_residues)

# -- FITTING-RELATED FUNCTIONS -----------

kf_means = np.array([-3.0] * num_residues)
kr_means = np.array([-3.0] * num_residues)
means = list(chain.from_iterable(zip(kf_means, kr_means)))
variances = np.array([1.0] * num_residues * 2)

def prior(mcmc, position):
    scaling_parameters = position[0:num_residues]
    scaling_prior = nbd_model_shared.prior(scaling_parameters,
                                           num_residues=num_residues)
    # Lognormal prior with given means and variances (in log10)
    return scaling_prior + \
           np.sum((position[num_residues:] - means)**2 / (2 * variances))


def random_initial_values(num_sets=1):
    initial_values_list = []

    kf_lower_bound = -5
    kf_upper_bound = -1
    kr_lower_bound = -5
    kr_upper_bound = -1

    for i in range(0, num_sets):
        kfs = 10 ** np.random.uniform(low=kf_lower_bound, high=kf_upper_bound,
                                      size=num_residues)
        krs = 10 ** np.random.uniform(low=kr_lower_bound, high=kr_upper_bound,
                                      size=num_residues)
        initial_values_list.append(list(chain.from_iterable(zip(kfs, krs))))

    # Concatenate with the randomized initial values from nbd_model_shared
    return np.concatenate((nbd_model_shared.random_initial_values(
                                        num_sets=num_sets, num_residues=5),
                           initial_values_list), axis=1)


def linear_pathway_groups(num_phases):
    """ Partition the sites into n groups.
    Rule 1: from all sites in solvent to insertion of sites in first group.
    Rule 2: from all sites in first group inserted to this plus all sites in second group.
    Rule 3: from all sites in first and second groups inserted to this plus all sites in third group
    Rule 4: ...etc."""

    # Should return list of (n) lists of sites (strings) 
    group_list = partition(Bax.sites, n)

    # Iterate over site groups
    all_prev = {}
    for i, group in enumerate(group_list):
        # Make the rule
        # Make the left hand side dict {'c3':'s', 'c62':'m'...}
        # Iterate over all previous groups
        for prev_group in enumerate(group_list[0:i]):
            # For each group, iterate over sites
            lhs_sites = [site for site in prev_group]
