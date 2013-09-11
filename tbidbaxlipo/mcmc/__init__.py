"""
Code to streamline the process of fitting mechanistic models by MCMC.

The bulk of this module consists of the class :py:class:`MCMC`, a subclass of
`bayessb.MCMC` (from the bayessb_ module) that contains code that is common to
all `tbidbaxlipo` MCMC fitting procedures (e.g., pickling/unpickling, plotting
fits vs. data, etc.).

In addition, this class contains a helper function,
:py:func:`import_mcmc_groups`, which parses lists of MCMC chain filenames into
group representing different runs (with different random seeds) of the same
model.

For more information on implementing fitting procedures by subclassing
:py:class:`tbidbaxlipo.mcmc.MCMC`, see
:ref:`creating_new_mcmc_fitting_procedures`.

.. _bayessb: http://sorgerlab.github.com/bayessb
"""

import bayessb
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import os

class MCMC(bayessb.MCMC):
    """Superclass for creating MCMC fitting procedures.

    The constructor chain-calls the `bayessb.MCMC` constructor, passing
    it the MCMCOpts object in ``options``.

    Parameters
    ----------
    options : bayessb.MCMCOpts
        The options for the MCMC fitting procedure.
    builder : :py:class:`tbidbaxlipo.models.core.Builder`
        An instance of :py:class:`tbidbaxlipo.models.core.Builder`
        (or a subclass).
    """
    def __init__(self, options, builder):
        bayessb.MCMC.__init__(self, options)

        self.builder = builder

    # Pickling functions for this class
    def __getstate__(self):
        """Pickling function.

        Because function references are not pickled correctly, this function
        sets the references to the likelihood, prior, and step functions
        (contained in ``self.options``) to None prior to pickling.
        """
        self.options.likelihood_fn = None
        self.options.prior_fn = None
        self.options.step_fn = None
        mcmc_state = bayessb.MCMC.__getstate__(self)
        return mcmc_state

    def __setstate__(self, state):
        """Unpickling function.

        Restores the references to the likelihood, prior, and step functions in
        ``self.options``. The likelihood function is contained by
        ``self.likelihood``, and the prior function is obtained from the model
        builder instance as ``self.builder.prior``. By default, the step
        function used is the generic step function contained in this class,
        namely ``self.step``.
        """
        bayessb.MCMC.__setstate__(self, state)
        self.options.likelihood_fn = self.likelihood
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    def do_fit(self):
        """Runs MCMC on the given model.

        Initializes the MCMC (by calling ``self.initialize()``), prints out the
        starting parameter values, then runs the MCMC (by calling
        ``self.run()``)
        """
        self.initialize()

        # Print initial parameter values
        init_vals = zip([p.name for p in self.options.model.parameters],
                          self.cur_params(position=self.initial_position))
        init_vals_str = 'Initial values:\n'
        init_vals_str += '\n'.join(['%s: %g' % (init_vals[i][0],
                                                 init_vals[i][1])
                                     for i in range(0, len(init_vals))])
        print "------------------------"
        print init_vals_str
        print "------------------------"

        # Run it!
        self.run()

    def fit_plotting_function(self, position):
        """Gets the observable timecourse and plots it against the data.

        Useful for visualizing model output by examining the fit of the model
        at the maximum likelihood or maximum a posteriori parameters, etc.

        Parameters
        ----------
        position : numpy.array
            Vector of parameter values to use for simulation and plotting.

        Returns
        -------
        matplotlib.figure.Figure
            Figure instance containing the plot of the data and the simulation
            at the given position.
        """
        # Run the simulation at the given position
        timecourses = self.get_observable_timecourses(position)
        # Make the plot
        fig = Figure()
        ax = fig.gca()
        # Add the data to the plot
        self.plot_data(ax)

        # Add the simulations to the plot
        for obs_name, timecourse in timecourses.iteritems():
            ax.plot(timecourse[0], timecourse[1], label=obs_name)

        # Label the plot
        ax.set_xlabel('Time')
        ax.set_ylabel('Concentration')
        fontP = FontProperties() 
        fontP.set_size('small')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='upper left', prop=fontP, ncol=1, bbox_to_anchor=(1, 1),
             fancybox=True, shadow=True)
        #ax.legend(loc='upper center', prop=fontP, ncol=5,
        #            bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True)
        canvas = FigureCanvasAgg(fig)
        fig.set_canvas(canvas)
        return fig

    @staticmethod
    def step(mcmc):
        """The function to call at every MCMC iteration.

        Prints out a few progress indicators such as the acceptance rate,
        likelihood, posterior, etc.

        Parameters
        ----------
        mcmc : bayessb.MCMC
            The MCMC being run.
        """
        window = 200
        local_acc = np.sum(mcmc.accepts[(mcmc.iter - window):mcmc.iter]) / \
                    float(window)

        if mcmc.iter % 20 == 0:
            print 'iter=%-5d  sigma=%-.3f  T=%-.3f  loc_acc=%-.3f  ' \
                  'glob_acc=%-.3f  lkl=%g  prior=%g  post=%g' % \
                  (mcmc.iter, mcmc.sig_value, mcmc.T,
                   local_acc,
                   mcmc.acceptance/(mcmc.iter+1.), mcmc.accept_likelihood,
                   mcmc.accept_prior, mcmc.accept_posterior)

    # "Virtual" functions that must be implemented
    def plot_data(self, axis):
        raise NotImplementedError()

    def get_observable_timecourses(self, position):
        raise NotImplementedError()

    def get_likelihood_function(self):
        raise NotImplementedError()

    def get_basename(self):
        raise NotImplementedError()

# Chain handling helper function
# ==============================
def import_mcmc_groups(filenames):
    """Loads the chains into groups representing multiple runs of the same
    model.

    Assumes that the pickle filenames are structured as

        basename = '%s_%s_%s_%d_s%d.xxx' % (model, nbd_sites, nbd_observables,
                                        nsteps, random_seed)

    With the suffix ``.xxx`` separated by a dot and the seed coming last in
    the underscore-separated arguments.

    Parameters
    ----------
    filenames : list of strings
        List of strings representing the chain filenames to be sorted into
        groups, e.g., of the type returned by ``glob.glob()``.

    Returns
    -------
    dict of lists of MCMC filenames
        The keys in the dict are the filename prefixes that represent the
        arguments to the MCMC procedure (e.g., ``tard_c3_iBax_4000``. Each dict
        entry contains a list of MCMC filenames associated with those run
        conditions.
    """

    mcmc_groups = {}

    for filename in filenames:
        # Split off the suffix from the filename
        (prefix, suffix) = filename.rsplit('.', 1)
        # Separate the filename into the final argument identifying the
        # random seed, and everything that comes before it:
        (mcmc_args, seed) = prefix.rsplit('_', 1)
        mcmc_args = os.path.basename(mcmc_args)

        # Check to see if we've imported another chain of this type already.
        # If so, add the current chain to the list of chains for this group:
        if mcmc_args in mcmc_groups:
            mcmc_groups[mcmc_args].append(filename)
        # If not, create a new entry in the dict containing this MCMC
        else:
            mcmc_groups[mcmc_args] = [filename]

    return mcmc_groups

