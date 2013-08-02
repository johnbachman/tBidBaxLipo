Submitting SSA Jobs
===================

It is often useful to run many stochastic simulations of pore formation so
that they can be aggregated and compared to a deterministic model. Doing so
generally requires that the jobs be run in parallel on a computing cluster. The
following is a protocol for how to do this using the scripts in this package.

**1. Write a script to build the SSA model and run it.**

The script should assemble the model programmatically using one
of the builder classes, e.g. :py:class:`tbidbaxlipo.models.n_cpt.Builder`.
Note that any parameter overrides to be specified via the ``params_dict``
argument should be specified here, as well as the scaling factor to use
for scaling agent numbers. For example::

    from tbidbaxlipo.models.n_cpt_jobs import Job
    from tbidbaxlipo.models.n_cpt import Builder

    scaling_factor = 10
    params_dict = {'Bax_transloc_kf':1e-2}

    b = Builder(scaling_factor=scaling_factor, params_dict=params_dict)
    b.build_model_bax_heat()

.. important::

    If you are planning on comparing the stochastic simulation output to a
    deterministic analog, make sure that you are setting up the stochastic and
    deterministic models in exactly the same way, including using the same set
    of parameters. Review any model-building motifs that have specific
    overrides for equivalence. In particular, `check the default values of any
    parameters that have not been overridden in the params_dict!`

In the script, create an instance of
:py:class:`tbidbaxlipo.models.n_cpt_jobs.Job`, setting the time range, the number
of simulation steps to record, and the number of simulations to run in a single
job::

    tmax = 10000
    n_steps = 200
    num_sims = 10
    j = Job(b.model, tmax, n_steps, num_sims)
    j.run()

You now have a script that will handle the execution of a single job (though
note that one job can run multiple simulations).

.. note:: Local execution

    Since there is nothing cluster-specific in the run script, you can also use
    it to do a test run locally.

**2. Submit the run script jobs via LSF.**

Log into the computing cluster and navigate to the directory on the server
(e.g. `tbidbaxlipo/simdata`) where you want the simulation results to be saved.
It is a good practice to put the run script used to create and run the model in
this directory so you maintain a link between the simulation setup and the
results. Use the script :py:mod:`tbidbaxlipo.models.n_cpt_jobs` to submit
the simulation jobs on Orchestra, as follows::

    python -m tbibaxlipo.models.n_cpt_jobs submit [run_script.py] [num_jobs]

**3. Parse the results back.**

In most cases you will only want the means and the standard deviations of the
observables at each time point. To calculate and save these, you run the script::

    python -m tbidbaxlipo.models.n_cpt_jobs parse /path/*.gdat

The `.gdat` files to load and parse are passed in as a glob or list of files.
The script saves the means of each observable in a record array pickled to
the file ``means.pck``; the standard deviations are pickled in ``stds.pck``.

**4. (Optional) Make the data an importable resource.**

If you put the simulation data in a submodule directory, you can add an
``__init__.py`` file which allows the mean and SD data to be imported for
plotting and analysis.

TBD

**4. Plot results and compare with deterministic model.**

TBD

