.. _creating_new_mcmc_fitting_procedures:

Creating new MCMC fitting procedures
====================================

A key aim of the code in the `tbidbaxlipo` module is to allow mechanistic
models of Bcl-2 interactions to be systematically calibrated to datasets by
Markov Chain Monte Carlo (MCMC) methods. The core MCMC algorithm used is
implemented in the (separate) bayessb_ module. However, there are a number of
steps that must be taken to adapt the MCMC code to work with particular types
of mechanistic models and particular types of experimental datasets (e.g., NBD
fluorescence vs. dye release, or single timecourse vs. titration).

.. _bayessb: http://sorgerlab.github.com/bayessb/

The process of writing fitting code for `tbidbaxlipo` mechanistic models is
streamlined by :py:class:`tbidbaxlipo.mcmc.MCMC`, a subclass of `bayessb.MCMC`.
This subclass contains a few methods for pickling and unpickling MCMC objects,
assigning the likelihood function and prior appropriately, and plotting model
fits against the data. Specific fitting procedures are then created by further
subclassing :py:class:`tbidbaxlipo.mcmc.MCMC` and implementing the appropriate
methods.

Examples of fitting procedures created in this way include:

* :py:class:`tbidbaxlipo.mcmc.nbd_mcmc.NBD_MCMC`, which is used for fitting
  NBD-Bax fluorescence data from the PTI fluorimeter;
* :py:class:`tbidbaxlipo.mcmc.nbd_multiconf.NBDPlateMCMC`, for fitting
  NBD-Bax fluorescence data from the plate reader;
* :py:class:`tbidbaxlipo.mcmc.pore_mcmc.PoreMCMC`, for fitting dye release
  titration data.

In addition to implementing the appropriate methods in the subclass, these
modules also contain boilerplate code for parsing and validating arguments from
the command line specifying parameters of the MCMC fitting process.

What follows is a step-by-step description of the process of creating a fitting
procedure for a particular combination of models/data.

1. Subclassing :py:class:`tbidbaxlipo.mcmc.MCMC`
------------------------------------------------

**Implement the __init__ method.**

After coming up with a good module and class name for your fitting procedure,
create the class, subclassing :py:class:`tbidbaxlipo.mcmc.MCMC`. To write the
``__init__`` method, you'll need to think about what kind of information your
fitting procedure will need to have access to. At a minimum the ``__init__``
method will need to take an instances of `bayessb.MCMCOpts` and
:py:class:`tbidbaxlipo.models.core.Builder` as arguments, so that these can be
passed to the ``__init__`` method of :py:class:`tbidbaxlipo.mcmc.MCMC`.

It will also be necessary to pass in the dataset being used for fitting, as
this will be required by the likelihood function, plotting function, etc.
Including a separate field for the name of the dataset is useful for creating
unique output filenames.

Finally, the ``__init__`` method should set the references to the likelihood,
prior, and step functions, as these are defined from the builder. Here
is an example from :py:class:`tbidbaxlipo.mcmc.pore_mcmc.PoreMCMC`::

    class PoreMCMC(tbidbaxlipo.mcmc.MCMC):
        def __init__(self, options, data, dataset_name, builder):
            # Call the superclass constructor
            tbidbaxlipo.mcmc.MCMC.__init__(self, options, builder)

            # Store the data
            self.data = data
            self.dataset_name = dataset_name

            # Set the MCMC functions
            self.options.likelihood_fn = self.likelihood
            self.options.prior_fn = self.builder.prior
            self.options.step_fn = self.step

**Implement the likelihood function.**

The likelihood function is the centerpiece of the MCMC script, where we define
the error metric (and hence the likelihood) for the model/data pair. The
function must ultimately be a function that takes two arguments, the MCMC
object (``mcmc``) and the position in parameter space (``position``).  It
should return a single value representing the likelihood of observing the data
given the model with the current set of parameters.

The likelihood function can be implemented in a variety of ways, including as a
closure around a function that is configured using information in the MCMC
object and then returned (as in
:py:meth:`tbidbaxlipo.mcmc.nbd_mcmc.NBD_MCMC.get_likelihood_function`) or as a
static method. Since the likelihood function is always called with two
arguments, the first of which is the MCMC object itself, any information in the
MCMC object will be accessible through this argument.

Likelihood functions involving a single simulation of the model are fairly
simple, involving running the model, calculated the chi-squared error relative
to the data, and returning the value. Datasets involving titrations are
slightly more complicated in that they involve multiple simulations of the model
with different initial conditions, possibly with different time vectors
(because each concentration condition in the dataset may have a different set
of time vectors).

The following is the implementation of the likelihood function for `PoreMCMC`,
:py:meth:`tbidbaxlipo.mcmc.pore_mcmc.PoreMCMC.likelihood`::

    @staticmethod
    def likelihood(mcmc, position):
        err = 0
        for bax_conc in mcmc.data.columns:
            # Get the data for this concentration
            tc = mcmc.data[bax_conc]
            y_data  = np.array(tc[:,'MEAN'])
            time = np.array(tc[:,'TIME'])
            mcmc.solver.tspan = time # set the time span

            # Get the simulated data for this concentration
            mcmc.options.model.parameters['Bax_0'].value = bax_conc
            x = mcmc.simulate(position=position, observables=True)
            avg_pores = x['pores']/ \
                        mcmc.options.model.parameters['Vesicles_0'].value
            y_mod = 1 - np.exp(-avg_pores)

            # Calculate the error, accounting for the SD at this
            # concentration.
            # Skip the first timepoint--the SD is 0 (due to normalization)
            # and hence gives nan when calculating the error.
            err += np.sum(((y_data[1:] - y_mod[1:])**2) / \
                    (2 * (np.array(tc[:,'SD'][1:]) ** 2)))
        return err

Note that the solver object contained by the MCMC instance must have its time
vector ``tspan``, and the initial condition for Bax ``Bax_0``, reset for each
simulation. The error at each concentration is calculated and the total error
is returned by the function.

**Implement the plot_data method.**

This function will be called by the superclass method
:py:meth:`tbidbaxlipo.mcmc.MCMC.fit_plotting_function`. It is used to plot the
data into a figure for comparing the fit of the model run with a set (or
multiple sets) of parameters. The data-plotting function for
:py:class:`tbidbaxlipo.mcmc.pore_mcmc` involves plotting each timecourse in the
titration, stored in a `pandas.Dataframe`::

    def plot_data(self, axis):
        # Plot the titration of Bax timecourses
        for bax_conc in self.data.columns:
            tc = self.data[bax_conc]
            axis.errorbar(tc[:,'TIME'], tc[:,'MEAN'], yerr=tc[:,'SD'],
                       color='gray')

**Implement the get_observable_timecourses method.**

This function takes a parameter vector and returns a dict containing the
simulated timecourses for the observables. The purpose of this function is
mainly for plotting model fits against the data--for example, it is called by
the superclass function
:py:meth:`tbidbaxlipo.mcmc.MCMC.fit_plotting_function`.

The form of this function will look something like that of the likelihood
function, as it will involving gathering up the results from the observables of
interest, possibly by iterating over a set of initial concentrations. The
results are returned in a somewhat unusual format: a dict of lists, where the
keys are the human-readable names for the observables or simulation conditions
(to be used in the plot legend), and the values are two-element lists
consisting of the time vector and the simulated values: ``[time, y]``.

As an example, here is the implementation for `PoreMCMC`, :py:meth:`tbidbaxlipo.mcmc.pore_mcmc.PoreMCMC.get_observable_timecourses`::

    def get_observable_timecourses(self, position):
        """Return the timecourses for all concentrations."""
        timecourses = collections.OrderedDict()

        for bax_conc in self.data.columns:
            # Get the timepoints for this concentration
            tc = self.data[bax_conc]
            time = np.array(tc[:,'TIME'])
            self.solver.tspan = time # set the time span
            self.options.model.parameters['Bax_0'].value = bax_conc
            x = self.simulate(position=position, observables=True)
            avg_pores = x['pores'] / \
                        self.options.model.parameters['Vesicles_0'].value
            y_mod = 1 - np.exp(-avg_pores)
            timecourses['Bax %d nM' % bax_conc] = [time, y_mod]
        return timecourses

**Implement the get_basename method.**

Finally, implement the ``get_basename`` method, which returns the string name
that will be used for pickled MCMC output files. The method should include
whatever information from the MCMC object that is necessary for the name to be
unique, such as the dataset and model used, the number of steps in the walk,
the random seed, etc. Here is the implementation for `PoreMCMC`::

    def get_basename(self):
        return '%s_%s_%s_%d_s%d' % (self.dataset_name,
                                 self.builder.get_module(),
                                 self.options.model.name,
                                 self.options.nsteps,
                                 self.options.seed)

2. Creating a run script
------------------------

After implementing the key fitting and plotting methods in the subclass of
:py:class:`tbidbaxlipo.mcmc.MCMC`, it is necessary to include a run script that
gets the dataset, parses arguments at the command-line, runs the MCMC, and
pickles the results. The run script can be implemented either in a separate
file or in a ``if __name__ == `__main__`:`` section in the top-level of the
module containing the rest of the code (the latter is the recommended
approach).

The code is for running the scripts is substantially boilerplate, but enough
differences exist regarding arguments for models and datasets to use that
generalizing the run script does not seem worthwhile. Instead, copy-and-paste
with modifications seems to be a satisfactory approach. Run script code can
be duplicated from :py:mod:`tbidbaxlipo.mcmc.pore_mcmc`, :py:mod:`tbidbaxlipo.mcmc.nbd_multiconf`, or :py:mod:`tbidbaxlipo.mcmc.nbd_mcmc_run`.

The run script, when implemented, should allow the execution of a MCMC
fitting procedure at the command-line using a syntax such as the following
(example from ``pore_mcmc.py``)::

    python -m tbidbaxlipo.mcmc.pore_mcmc random_seed=0 model=bax_heat \
                cpt_type=one_cpt nsteps=1000

3. Creating a submission script
-------------------------------

The purpose of the submission script is to streamline the process of submitting
many parallel MCMC jobs on the Orchestra computing cluster. Because the types
of models, data, temperatures, or other parameters to systematically iterate
over for job submission may vary depending on the type of model, each
model/data type will likely require its own job submission script. However, the
implementation tends to follow a fairly boilerplate pattern, with the sets to
iterate over defined as a number of lists which are then iterated over; each
combination of parameters is then executed as a distinct job by a call to a
command-line operation.

Rather than a large set of nested for loops, an improved approach is to take
the Cartesian product of the various lists using `itertools.product`, and then
iterating over the result. Here is an example submission loop from
:py:mod:`tbidbaxlipo.mcmc.pore_mcmc_jobs`, the job submission script
corresponding to `pore_mcmc`::

    # Iterate over the Cartesian product of the different argument lists
    for args in itertools.product(model_arg_list,
                                  cpt_type_arg_list,
                                  random_seed_arg_list):
        fixed_args = ['nsteps=%d' % nsteps]
        cmd_list = base_cmd_list(output_filename_from_args(args)) + \
                   list(args) + fixed_args
        print ' '.join(cmd_list)
        subprocess.call(cmd_list)

4. Creating a parallel tempering script
---------------------------------------

