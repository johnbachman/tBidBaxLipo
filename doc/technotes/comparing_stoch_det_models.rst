Comparing stochastic and deterministic simulations
==================================================

It is often useful to run stochastic simulations of pore formation so that they
can be aggregated and compared to a deterministic model. Doing so requires that
the stochastic and deterministic models be built up in precisely the same
fashion, with the same topology and parameters. Running a sufficient number of
stochastic simulations to perform a meaningful comparison also often requires
that the simulations be run in parallel on a computing cluster. The following
is a protocol for how to do this using the scripts in this package.

1. Write a script to build the SSA model and run it.
----------------------------------------------------

The module :py:mod:`tbidbaxlipo.models.simulation` contains most of the
boilerplate code involved in setting up and parsing stochastic simulation
results. To create a script for simulating a specific model, create a new file
that will be your run script. In this script you will create a subclass of
:py:class:`tbidbaxlipo.models.simulation.Job`, implementing only the methods
``__init__`` and ``build``. The script should also declare top-level package
variables ``jobs`` and ``job_name``. The list variable ``jobs`` will contain
the instances of the ``Job`` class (initialized with various initial
conditions, parameters, etc.) that are to be run. Finally you will add a
boilerplate ``__main__`` script for running the job at the command-line.  Here
is a complete example::

    from tbidbaxlipo.models import simulation

    class Job(simulation.Job):
        def __init__(self, Bax_transloc_kf):
            params_dict = {'Bax_transloc_kf':Bax_transloc_kf}
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

    jobs = [Job(1e-2), Job(1e-3)]

    if __name__ == '__main__':
        simulation.run_main(jobs)

The ``__init__`` method's only role is to make explicit the parameter
overrides, duration, and other features of the simulation that are passed
to the superclass constructor.

The ``build`` method takes a module (e.g. :py:mod:`tbidbaxlipo.models.one_cpt`,
:py:mod:`tbidbaxlipo.models.n_cpt` or :py:mod:`tbidbaxlipo.models.site_cpt`)
and returns the appropriate builder instance after the model has been
constructed.  The ``build`` method should call the constructor
``module.Builder`` as shown, passing in the ``params_dict`` and
``scaling_factor`` arguments (note that ``scaling_factor`` will simply be
ignored by the ``__init__`` method of ``one_cpt.Builder``).  Building the model
this way requires that model assembly code is written in one place and the same
parameters are passed in.

.. important::

    For comparing deterministic and stochastic versions of models, make sure
    that any module-specific implementations of model building motifs are
    equivalent. In particular, `check the default values of any parameters that
    have not been overridden in the params_dict!`

Finally, if stochastic simulations are to be run on the computing cluster, it
is necessary to add the final ``__main__`` script, which calls the generic
``run_main`` method from the ``simulation`` module. ``run_main`` handles
parsing of command line arguments and running the appropriate job in the
``jobs`` list based on the command line arguments.

.. note:: Local execution

    Since there is nothing cluster-specific in the run script, you can also use
    it to do a test run locally. However, if there is more than one job
    in the job list, you will have to pass in an index.

2. Submit the run script jobs via LSF.
--------------------------------------

Log into the computing cluster and navigate to the directory on the server
(e.g. `tbidbaxlipo/plots/mysim`) where you want the simulation results to be
saved. It is a good practice to write the run script used to create and run the
model as the ``__init__.py`` file for the module ``mysim``. This maintains the
link between the simulation setup and the results. Run the script
:py:mod:`tbidbaxlipo.models.simulation` at the command-line to submit the
simulation jobs on Orchestra, as follows::

    python -m tbidbaxlipo.models.simulation submit [run_script.py] [num_jobs]

The submit script will submit ``num_jobs`` jobs to Orchestra for `each` of the
jobs in the job array.

Jobs run in this way will automatically save simulation results in a series of
directories ``data_0``, ``data_1``, etc. The numbers associated with the
directories are matched to the index of the jobs in the ``jobs`` list of the
run script.

3. Parse the results back.
--------------------------

The parsing functionality of :py:mod:`tbidbaxlipo.models.simulation` will parse
the simulation results into an HDF5 file. There are two general approaches to
parsing: the first is to load all of the BNG files and write the final HDF5
file in one process; the second is to pre-parse BNG files into assemble
separate HDF5 files for each simulation condition (each Job instance in the job
list) and then assemble them into the final HDF5 file. The latter approach has
the advantage that the rate-limiting step of BNG file parsing can be run in
parallel.

parse_single
~~~~~~~~~~~~

To parse data from a single job/condition, run::

    python -m tbidbaxlipo.models.simulation parse_single [hdf5_file] [dir]

For example::

    python -m tbidbaxlipo.models.simulation parse_single data data_0

This will parse the .gdat files in the directory data_0 into the file
data.hdf5.

parse_set
~~~~~~~~~

To parse data from a list of different conditions, run::

    python -m tbidbaxlipo.models.simulation parse_set [hdf5_file] [dir_base] [num_dirs]

This will iteratively parse .gdat files from the directories '%s%d' %
(dir_base, i) for i ranging from 0 through num_dirs - 1.

For example::

    python -m tbidbaxlipo.models.simulation parse_set data data_ 20

Note that if there are many simulations or many conditions, the parsing process
can take a long time.

parse_parallel
~~~~~~~~~~~~~~

Since parse_set can take a very long time to iterate over all of the .gdat
files, it is often more efficient to submit a parallel job to parse the data in
stages. To do this, open a shell on Orchestra and run::

    python -m tbidbaxlipo.models.simulation parse_parallel [hdf5_file] [dir_base] [num_dirs]

Note that the arguments are the same as for parse_set. This will submit an LSF
job to parse each of the data directories (using parse_single). The result
will be a set of HDF5 files with names corresponding to the data directories,
e.g. data_0.hdf5, data_1.hdf5, etc. These HDF5 files can then be assembled into
the full HDF5 file using parse_assemble.

parse_assemble
~~~~~~~~~~~~~~

Run this after generating a set of HDF5 files via parse_parallel. The syntax
is similar to parse_set, except the arguments refer to HDF5 files rather
than data directories::

    python -m tbidbaxlipo.models.simulation parse_assemble [hdf5_file] [hdf5_base] [num_files]

For example, to parse 20 HDF5 files data_0.hdf5, data_1.hdf5, etc. into a single HDF5 output file data, run::

    python -m tbidbaxlipo.models.simulation parse_assemble data data_ 20

Note here that the .hdf5 suffix is added automatically.

4. (Optional) Make the data an importable resource.
---------------------------------------------------

If you put the simulation data in a submodule directory, you can add a few
boilerplate lines to the ``__init__.py`` file which will allow access to the
data in the HDF5 file::

    mod_path = os.path.dirname(sys.modules[__name__].__file__)
    hdf5_filename = os.path.abspath(os.path.join(mod_path, 'data.hdf5'))
    if os.path.exists(hdf5_filename):
        data = simulation.CptDataset(hdf5_filename)

If the name ``data.hdf5`` has been used for the HDF5 file, then this code
can be used without modification.

The HDF5 dataset can then be imported by calling code as::

    from tbidbaxlipo.plots.mysim import data

5. Plot results and compare with deterministic model.
-----------------------------------------------------

.. todo:: The below needs to be updated.

The class :py:class:`tbidbaxlipo.models.simulation.Job` contains a
:py:meth:`tbidbaxlipo.models.simulation.Job.run_one_cpt` that handles the
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
