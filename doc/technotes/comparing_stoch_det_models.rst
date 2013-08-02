Comparing stochastic and deterministic simulations
==================================================

It is often useful to run stochastic simulations of pore formation so that they
can be aggregated and compared to a deterministic model. Doing so requires that
the stochastic and deterministic models be built up in precisely the same
fashion, with the same topology and parameters. Running a sufficient number of
stochastic simulations to perform a meaningful comparison also often requires
that the simulations be run in parallel on a computing cluster. The following
is a protocol for how to do this using the scripts in this package.

**1. Write a script to build the SSA model and run it.**

The module :py:mod:`tbidbaxlipo.models.n_cpt_jobs` contains most of the
boilerplate code involved in setting up and parsing stochastic simulation
results. To create a script for simulating a specific model, create a new file
that will be your run script. In this script you will create a subclass of
:py:class:`tbidbaxlipo.models.n_cpt_jobs.Job`, implementing only the methods
``__init__`` and ``build``. In addition you will add a short main script for
running the job at the command-line.  Here is a complete example::

    from tbidbaxlipo.models import n_cpt_jobs

    class Job(n_cpt_jobs.Job):
        def __init__(self):
            params_dict = {'Bax_transloc_kf':1e-2}
            scaling_factor = 5
            tmax = 10000
            n_steps = 200
            num_sims = 3
            super(Job, self).__init__(params_dict, scaling_factor, tmax,
                                      n_steps, num_sims)

        def build(self, module):
            builder = module.Builder(params_dict=self.params_dict,
                                     scaling_factor=self.scaling_factor)
            builder.build_model_bax_heat()
            return builder

    if __name__ == '__main__':
        j = Job()
        j.run_n_cpt()

The ``__init__`` method's only role is to make explicit the parameter
overrides, duration, and other features of the simulation that are passed
to the superclass constructor.

The ``build`` method takes a module (generally either
:py:mod:`tbidbaxlipo.models.n_cpt` or :py:mod:`tbidbaxlipo.models.one_cpt` and
returns the appropriate builder instance after the model has been constructed.
The ``build`` method should call the constructor ``module.Builder`` as shown,
passing in the ``params_dict`` and ``scaling_factor`` arguments (note that
``scaling_factor`` will simply be ignored by the ``__init__`` method of
``one_cpt.Builder``). Building the model this way requires that model assembly
code is written in one place and the same parameters are passed in.

.. important::

    For comparing deterministic and stochastic versions of models, make sure
    that any module-specific implementations of model building motifs are
    equivalent. In particular, `check the default values of any parameters that
    have not been overridden in the params_dict!`

Finally, if stochastic simulations are to be run on the computing cluster, it
is necessary to add the final ``__main__`` script, which simply creates an
instance of the ``Job`` class and calls the ``run_n_cpt`` method. You now have
a script that will handle the execution of a single job (though note that one
job can run multiple simulations).

.. note:: Local execution

    Since there is nothing cluster-specific in the run script, you can also use
    it to do a test run locally.

**2. Submit the run script jobs via LSF.**

Log into the computing cluster and navigate to the directory on the server
(e.g. `tbidbaxlipo/simdata/mysim`) where you want the simulation results to be
saved. It is a good practice to put the run script used to create and run the
model in this directory so you maintain a link between the simulation setup and
the results. Run the script :py:mod:`tbidbaxlipo.models.n_cpt_jobs` at the
command-line to submit the simulation jobs on Orchestra, as follows::

    python -m tbidbaxlipo.models.n_cpt_jobs submit [run_script.py] [num_jobs]

**3. Parse the results back.**

In most cases you will only want the means and the standard deviations of the
observables at each time point. To calculate and save these, you run the
script::

    python -m tbidbaxlipo.models.n_cpt_jobs parse /path/*.gdat

The `.gdat` files to load and parse are passed in as a glob or list of files.
The script saves the means of each observable in a record array pickled to
the file ``means.pck``; the standard deviations are pickled in ``stds.pck``.

**4. (Optional) Make the data an importable resource.**

If you put the simulation data in a submodule directory, you can add an
``__init__.py`` file which allows the mean and SD data to be imported for
plotting and analysis. Here is some boilerplate code::

    import pickle
    import pkgutil

    try:
        means = pickle.loads(pkgutil.get_data(
                            'tbidbaxlipo.simdata.sim_test', 'means.pck'))
        stds = pickle.loads(pkgutil.get_data(
                            'tbidbaxlipo.simdata.sim_test', 'stds.pck'))
    except IOError:
        pass

The try/catch block handles the case when the pickle files don't exist (and
hence allows submodules to be imported without error).  Only the path to the
package containing the data (``tbidbaxlipo.simdata.sim_test``) needs to
be changed from the above example to use it for a new dataset.

**5. Plot results and compare with deterministic model.**

The class :py:class:`tbidbaxlipo.models.n_cpt_jobs.Job` contains a
:py:meth:`tbidbaxlipo.models.n_cpt_jobs.Job.run_one_cpt` that handles the
construction and simulation of ``one_cpt`` models in precisely analogous
fashion to ``n_cpt``, streamlining comparison of models. Here is an example
plotting script::

    from matplotlib import pyplot as plt
    from tbidbaxlipo.simdata.sim_test import means, stds
    from tbidbaxlipo.simdata.sim_test.run_script import Job

    # Create the job instance
    j = Job()

    # Run the deterministic simulation
    (t, det_obs) = j.run_one_cpt()

    # Plot deterministic results
    plt.ion()
    plt.figure()
    plt.plot(t, det_obs['pores'])

    # Plot stochastic results
    plt.errorbar(means['time'], means['pores'] / j.scaling_factor,
                 yerr=stds['pores'] / j.scaling_factor)

In this example, note:

- ``means`` and ``stds`` are imported by using the resource strategy described
  above.
- An instance of ``run_script.Job`` is created to get access to the
  ``run_one_cpt`` method for deterministic simulation.
- If a scaling factor was used for stochastic simulation, rescaling of the
  observables may be required. Here, the instance of ``Job`` contains the
  scaling factor that was used for stochastic simulation, and hence it can
  be used to rescale the observables in ``means`` and ``stds``.

Running model comparisons locally
---------------------------------

In some cases SSA simulations run relatively fast and can be run locally rather
than on the computing cluster. This makes steps **2** - **4** above
unnecessary.  Instead, SSA results can be obtained and plotted within a single
script.  The following is an analogous example to **5** above::

    from matplotlib import pyplot as plt
    from tbidbaxlipo.simdata.sim_test.run_script import Job

    # Create the job instance
    j = Job()

    # Run the deterministic simulation
    (t, det_obs) = j.run_one_cpt()

    # Plot deterministic results
    plt.ion()
    plt.figure()
    plt.plot(t, det_obs['pores'], label='one_cpt')

    # Run the stochastic simulation
    xrecs = j.run_n_cpt(cleanup=True)
    (means, stds) = j.calculate_mean_and_std(xrecs)

    # Plot stochastic results
    plt.errorbar(means['time'], means['pores'] / j.scaling_factor,
                 yerr=stds['pores'] / j.scaling_factor, label='n_cpt')

    # Label the plot
    plt.xlabel('Time (secs)')
    plt.ylabel('Total pores')
    plt.title('Comparing one_cpt and n_cpt simulations')
    plt.legend(loc='lower right')

When run with the example run script shown above, this script produces the
following results:

.. plot::

    import tbidbaxlipo.simdata.sim_test.plot_comparison_local

