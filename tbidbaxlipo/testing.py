import numpy as np
from bayessb.convergence import convergence_criterion

class NBD_Tests(object):
    """.. todo:: document this class."""

    def __init__(self, mcmc, mask=None, thin=1):
        """.. todo:: document this function."""
        self.mcmc = mcmc
        self.results = []

        if mask == None:
            mask = self.mcmc.options.nsteps / 2

        (mixed_accepts, mixed_accept_steps) = mcmc.get_mixed_accepts(mask, thin)
        self.mixed_accepts = mcmc.get_mixed_accepts
        self.mixed_accept_steps = mixed_accept_steps

        self.diagnostic_results = []
        self.descriptive_results = []
        self.constraint_results = []

        self.example_timecourse = self.mcmc.simulate(position=mixed_accepts[0],
                                                     observables=True)
    
    # Diagnostic tests
    # ================
    def maximum_likelihood(self):
        return np.min(self.mcmc.likelihoods[self.mixed_accept_steps])

    def maximum_posterior(self):
        return np.min(self.mcmc.posteriors[self.mixed_accept_steps])

    def parameter_convergence_range(self):
        """.. todo:: Need a way to adapt this to several chains."""
        pass # TODO

    def run_diagnostic_tests(self):
        #import ipdb; ipdb.set_trace()
        self.diagnostic_results.append(('Maximum likelihood',
                                        self.maximum_likelihood())) 
        self.diagnostic_results.append(('Maximum posterior',
                                        self.maximum_posterior()))

    # Descriptive tests
    # =================
    def has_iBax(self):
        """.. todo:: better way to do all of these would be to check observables
        against list of species."""
        return np.any(self.example_timecourse['iBax'] != 0)

    def has_aBax(self):
        """.. todo:: better way to do all of these would be to check observables
        against list of species."""
        return np.any(self.example_timecourse['aBax'] != 0)

    def has_Bax2(self):
        """.. todo:: better way to do all of these would be to check observables
        against list of species."""
        return np.any(self.example_timecourse['Bax2'] != 0)

    def run_descriptive_tests(self):
        self.descriptive_results.append(('Has iBax?',
                                         self.has_iBax()))
        self.descriptive_results.append(('Has aBax?',
                                         self.has_aBax()))
        self.descriptive_results.append(('Has Bax2?',
                                         self.has_Bax2()))

    # Phenomenological tests/constraints
    # ==================================
    def tBid_Bax_monotonically_increasing(self):
        """ .. todo:: document the basis for this"""
        return monotonic_increasing(self.example_timecourse['tBidBax'])

    def iBax_monotonically_increasing(self):
        """ .. todo:: document the basis for this"""
        return monotonic_increasing(self.example_timecourse['iBax'])

    def run_constraints(self):
        """.. todo:: document concept of constraint tests"""
        self.constraint_results.append(('tBid/Bax monotonically increasing?',
                                    self.tBid_Bax_monotonically_increasing()))
        self.constraint_results.append(('iBax monotonically increasing?',
                                    self.iBax_monotonically_increasing()))

# Helper functions
# ================
def monotonic_increasing(x):
    # TODO rewrite so doesn't allow fixed, unchanging values
    dx = np.diff(x)
    return np.all(dx >= 0)
           
def monotonic_decreasing(x):
    # TODO rewrite so doesn't allow fixed, unchanging values
    dx = np.diff(x)
    return np.all(dx <= 0)
        
