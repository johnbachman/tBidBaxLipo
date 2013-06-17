"""
A class for fitting the mechanistic pore formation models to dye release data.
"""

class PoreMCMC(tbidbaxlipo.mcmc.MCMC):
    """Fit mechanistic tBid/Bax models to dye release titration data.

    Initialize parent tbidbaxlipo.mcmc.MCMC and then set additional fields.

    Parameters
    ----------
    options : MCMCOpts
        Options for MCMC initialization.
    data :
        pandas data structure containing the timecourses for each Bax condition.
    builder : tbidbaxlipo.models.core.Builder
    """
    def __init__(self, options, data, builder):
        # Call the superclass constructor
        tbidbaxlipo.mcmc.MCMC.__init__(self, options, builder)

        # Set the MCMC functions
        self.options.likelihood_fn = self.get_likelihood_function()
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    # Implementations of necessary functions
    def get_likelihood_function(self):
        return None

    def plot_data(self, axis):
        raise NotImplementedError()

    def get_observable_timecourses(self, position):
        raise NotImplementedError()

    def get_basename(self):
        return '%s_%s_%d_s%d' % (self.dataset_name,
                                 self.options.model.name,
                                 self.options.nsteps,
                                 self.options.seed)


