from pysb import bng, kappa
import sys
import numpy as np
import pickle
from tbidbaxlipo.models import one_cpt, n_cpt, site_cpt
from pysb.integrate import Solver
import collections
import subprocess
import h5py
import os
import glob

TIME_CHUNK_SIZE = 101
"""The chunk size for the time dimension of saved HDF5 datafiles."""

SIM_CHUNK_SIZE = 100
"""The chunk size for the simulation dimension of saved HDF5 datafiles."""

class Job(object):
    """Class for storing aspects of a SSA simulation job.

    Parameters
    ----------
    params_dict : dict
        dict mapping parameter names to values for use in overriding
        default values during model construction.
    scaling_factor : int
        Number used for scaling the number of agents used in stochastic
        simulations. Ignored in deterministic model construction.
    tmax : number
        The duration of the simulation.
    n_steps : int
        The number of simulation data points to save and return.
    num_sims : int
        The number of simulations to run as part of this job.
    """
    def __init__(self, params_dict, scaling_factor, tmax, n_steps,
                 num_sims):
        self.params_dict = params_dict
        self.scaling_factor = scaling_factor
        self.tmax = tmax
        self.n_steps = n_steps
        self.num_sims = num_sims
        self._one_cpt_builder = None
        self._n_cpt_builder = None
        self._site_cpt_builder = None

    def build(self, module):
        """Virtual method to build the model from the appropriate model Builder.

        Must be implemented by any subclass. Should take a module as object
        (e.g., :py:mod:`tbidbaxlipo.models.one_cpt` or
        :py:mod:`tbidbaxlipo.models.n_cpt`) and return an instance of a Builder
        class after model construction is completed.  """

        raise NotImplementedError()

    def one_cpt_builder(self):
        """Lazy initialization method for the one_cpt builder."""
        if self._one_cpt_builder is None:
            self._one_cpt_builder = self.build(one_cpt)
        return self._one_cpt_builder

    def n_cpt_builder(self):
        """Lazy initialization method for the n_cpt builder."""
        if self._n_cpt_builder is None:
            self._n_cpt_builder = self.build(n_cpt)
        return self._n_cpt_builder

    def site_cpt_builder(self):
        """Lazy initialization method for the site_cpt builder."""
        if self._site_cpt_builder is None:
            self._site_cpt_builder = self.build(site_cpt)
        return self._site_cpt_builder

    def run_n_cpt(self, cleanup=False):
        """Run a set of simulations using the n_cpt implementation.

        Builds the model using the :py:meth:`Job.build` method, then runs
        the number of simulations specified in ``self.num_sims`` using
        `pysb.bng.run_ssa` and returns the results.

        Parameters
        ----------
        cleanup : boolean
            If True, specifies that the files created by BNG during simulation
            should be deleted after completion. Defaults is False.

        Returns
        -------
        list of numpy.recarrays
            List of record arrays, each one containing the results of a
            single stochastic simulation. The entries in the record array
            correspond to observables in the model.
        """
        b = self.n_cpt_builder()
        xrecs = []
        for i in range(self.num_sims):
            print "Running BNG simulation %d of %d..." % (i+1, self.num_sims)
            xrecs.append(bng.run_ssa(b.model, t_end=self.tmax,
                         n_steps=self.n_steps, cleanup=cleanup, output_dir='.',
                         output_file_basename='%s_proc%d_sim%d' %
                                    (b.model.name, os.getpid(), i)))
        return xrecs

    def run_site_cpt(self):
        """Run a set of simulations using the site_cpt implementation.

        Builds the model using the :py:meth:`Job.build` method, then runs
        the number of simulations specified in ``self.num_sims`` using
        `pysb.kappa.run_simulation` and returns the results.

        Returns
        -------
        list of numpy.recarrays
            List of record arrays, each one containing the results of a
            single stochastic simulation. The entries in the record array
            correspond to observables in the model.
        """
        b = self.site_cpt_builder()
        xrecs = []
        for i in range(self.num_sims):
            xrecs.append(kappa.run_simulation(b.model, time=self.tmax,
                                              points=self.n_steps,
                                              output_dir='.'))
        return xrecs

    def run_one_cpt(self):
        """Run a deterministic simulation using the one_cpt implementation.

        Builds the model using the :py:meth:`Job.build` method, then runs
        a deterministic simulation for the duration specified by
        ``self.tmax`` and returns the results.

        Returns
        -------
        tuple of (numpy.array, numpy.recarray)
            The first entry is the vector of timepoints; the second is the
            record array of the model observables.
        """

        b = self.one_cpt_builder()
        # We add one extra point so that the time coordinates of the
        # one_cpt and n_cpt/site_cpt simulations line up
        t = np.linspace(0, self.tmax, self.n_steps + 1)
        s = Solver(b.model, t)
        s.run()
        return (t, s.yobs)

class CptDataset(object):
    """Wrapper for HDF5 datasets of multi-compartment simulations.

    Attributes
    ----------
    self.datafile : the h5py File object
    self.sim_data : h5py dataset, indexed by (cond, sim, obs, t)
    self.dtypes : dtypes
    self.obs_dict : dict
    """
    def __init__(self, filename, data_name='data', dtypes_name='dtypes'):
        self.datafile = h5py.File(filename)
        self.sim_data = self.datafile[data_name]
        self.dtypes = pickle.loads(self.datafile[dtypes_name][:])
        self.obs_dict = dict((name, i)
                             for i, name in enumerate(self.dtypes.names))
        self._means = None
        self._sds = None

    def means(self, cond_index, obs_name):
        """Gets the mean timecourse for the given condition and observable.

        The mean is calculated across the simulations on a per-condition basis,
        so each condition has its own set of timecourse means for all of its
        observables.

        Parameters
        ----------
        cond_index : int
            Index for the simulation condition for which we want the means.
        obs_name : string
            The name of the observable.

        Returns
        -------
        numpy.array
            One-dimensional array with the timecourse.
        """

        # Lazy initialization of the timecourse means
        if self._means is None:
            self._means = np.mean(self.sim_data, axis=1)
        return self._means[cond_index, self.obs_dict[obs_name], :]

    def sds(self, cond_index, obs_name):
        """Gets the SD of the timecourse for the given condition and observable.

        The standard deviationis calculated across the simulations on a
        per-condition basis, so each condition has its own set of timecourse
        SDs for all of its observables.

        Parameters
        ----------
        cond_index : int
            Index for the simulation condition for which we want the standard
            deviations.
        obs_name : string
            The name of the observable.

        Returns
        -------
        numpy.array
            One-dimensional array with the timecourse.
        """
        import ipdb; ipdb.set_trace()

        # Lazy initialization of the timecourse means
        if self._sds is None:
            self._sds = np.std(self.sim_data, axis=1)
        return self._sds[cond_index, self.obs_dict[obs_name], :]

    def means_iterative(self, cond_index, obs_name):
        """Same as :py:method:`means` but prevents loading the entire
        (uncompressed) dataset into memory by iterating over slices of the
        matrix.
        """
        (n_conds, n_reps, n_obs, n_tpts) = self.sim_data.shape
        if self._means is None:
            self._means = np.empty((n_conds, n_obs, n_tpts))
            for cond_ix in xrange(n_conds):
                for obs_ix in xrange(n_obs):
                    self._means[cond_ix, obs_ix, :] = \
                            np.mean(self.sim_data[cond_ix, :, obs_ix, :])
                print '%d' % cond_ix
        return self._means[cond_index, self.obs_dict[obs_name], :]

    def sds_iterative(self, cond_index, obs_name):
        """Same as :py:method:`sds` but prevents loading the entire
        (uncompressed) dataset into memory by iterating over slices of the
        matrix.
        """
        (n_conds, n_reps, n_obs, n_tpts) = self.sim_data.shape
        if self._sds is None:
            self._sds = np.empty((n_conds, n_obs, n_tpts))
            for cond_ix in xrange(n_conds):
                for obs_ix in xrange(n_obs):
                    self._sds[cond_ix, obs_ix, :] = \
                            np.std(self.sim_data[cond_ix, :, obs_ix, :])
                print '%d' % cond_ix
        return self._sds[cond_index, self.obs_dict[obs_name], :]

    def get_mean_dye_release(self, cond_index, pore_obs_prefix='pores'):
        """Get the mean timecourse for dye release across a set of simulations.

        Parameters
        ----------
        cond_index : int
            The index of the condition to get the dye release timecourse for.
        pore_obs_prefix : string
            The observables denoting pores in each compartment are expected to
            follow a naming scheme with a common prefix followed by the name of
            the compartment, etc., e.g. 'pores_c1', 'pores_c2', etc.  The
            ``pore_obs_prefix`` argument specifies the common string prefix to
            the compartment pore observables. Default is "pores".

        Returns
        -------
        tuple of numpy.arrays
            The first element in the tuple is a numpy.array containing the mean
            timecourse for dye release; the second element contains the
            standard deviation across for dye release across the set of
            simulations.
        """
        # Get the list of pore observables; makes the assumption that all pore
        # observable names have the same format, with a prefix followed by the
        # compartment identifier, e.g. 'pores_c38'.
        pore_indices = [i for i, name in enumerate(self.dtypes.names)
                       if name.startswith(pore_obs_prefix + '_')]
        # Slice the data to get just the condition and pore observables.
        # The slice has shape (num_simulations, num_cpts, num_timepoints)
        data_slice = self.sim_data[cond_index, :, pore_indices, :]
        # Create a matrix for dye release with shape (num_sims, num_tpts)
        num_simulations = data_slice.shape[0]
        num_timepoints = data_slice.shape[2]
        num_vesicles = len(pore_indices)
        dye_release = np.zeros((num_simulations, num_timepoints))
        # For every simulation result, calculate a dye release vector
        for sim_index in range(num_simulations):
            for timepoint in range(num_timepoints):
                dye_release[sim_index, timepoint] = \
                      (np.count_nonzero(data_slice[sim_index, :, timepoint]) /
                       float(num_vesicles))
        # Now, take the average/SD across all of the simulations
        return (np.mean(dye_release, axis=0), np.std(dye_release, axis=0))

    def get_frequency_matrix(self, cond_index, obs_basename, timepoint):
        """Get the frequencies of observable values across the simulations.

        The primary use is to calculate expected distributions of molecules
        across compartments from many stochastic multi-compartment simulations.
        Here ``observable_names`` would be the list of observables for the
        molecule at the various compartments.

        Parameters
        ----------
        observable_names : list of string
            The list of observable names to look for in the simulation output
            array, e.g. 'pores_c1', 'pores_c2', etc.
        output : numpy record array
            A record array of simulation results indexed by observable names.
        timepoint : int
            The index of the timepoint to get the snapshot of the observable.
            If no value is provided, the final timepoint is used.

        Returns
        -------
        tuple of numpy.arrays
            The first entry is an array of indices giving the values over which
            the frequencies are observed. The second entry is a matrix with
            shape (len(index), len(output)), where (i, j) gives the number of
            times the the value index(i) was observed across the list of
            observables in simulation output(j)."""

        # Get the indices of the observables
        obs_indices = [i for i, name in enumerate(self.dtypes.names)
                       if name.startswith(obs_basename + '_')]
        # Slice the data to get just the condition, observables and timepoint
        # of interest
        data_slice = self.sim_data[cond_index, :, obs_indices, timepoint]
        # The list of counts
        counts_list = []
        for sim_index in range(data_slice.shape[0]):
            # Convert to a dict of frequencies of each amount
            values = data_slice[sim_index, :]
            counts_list.append(collections.Counter(values))

        # Get the set containing the amounts observed across all simulations
        all_keys = set()
        for counts in counts_list:
            all_keys |= set(counts.keys())

        # Convert to a matrix in which the 0th entry represents the frequency
        # of observing the smallest observed amount, and the last entry
        #represents the frequency of observing the largest observed amount
        key_min = min(all_keys)
        key_max = max(all_keys)
        freq_matrix = np.zeros((key_max - key_min + 1, len(counts_list)))
        for i, counts in enumerate(counts_list):
            for key, val in counts.iteritems():
                freq_matrix[key - key_min, i] = val
        # The index runs from the minimum to the maximum amount
        index = np.array(range(int(key_min), int(key_max) + 1))
        return (index, freq_matrix)

    def get_means_across_cpts(self, cond_index, obs_basename, timepoint):
        # Get the indices of the observables
        obs_indices = [i for i, name in enumerate(self.dtypes.names)
                       if name.startswith(obs_basename + '_')]
        # Slice the data to get just the condition, observables and timepoint
        # of interest
        data_slice = self.sim_data[cond_index, :, obs_indices, timepoint]

        obs_means = np.mean(data_slice, axis=1)
        obs_vars = np.var(data_slice, axis=1)
        return (obs_means, obs_vars)

def save_bng_dirs_to_hdf5(data_dirs, filename):
    """Load the .gdat files in each of the listed directories to HDF5.

    The HDF5 file created contains two datasets:

    * A dataset named ``'data'`` which contains the simulation results in a
      four-dimensional numpy array of floats with shape `(num_conditions,
      num_simulations, num_observables, num_timepoints)`.

    * A dataset named ``'dtypes'`` containing the string-pickled dtype of the
      simulation results. The dtype contains the names of all of the
      observables from the original record array, and hence can be used to
      get the index for a particular observable in the full table of simulation
      results.

    Parameters
    ----------
    data_dirs : list of strings
        The list of directories containing the .gdat files to load, e.g.
        `['data_0', 'data_1', 'data_2', ... ]`
    filename : string
        The name of the HDF5 file to create.

    Returns
    -------
    h5py dataset
        The dataset containing the parsed simulation results and the corresponding
        dtype list.
    """
    num_conditions = len(data_dirs)
    dataset = None

    print "Loading BNG simulation results from data directories..."
    for condition_index, data_dir in enumerate(data_dirs):
        print "Loading data from directory %s" % data_dir
        gdat_files = glob.glob('%s/*.gdat' % data_dir)

        for sim_index, gdat_file in enumerate(gdat_files):
            ssa_result = bng._parse_bng_outfile(gdat_file)
            # Initialize the dataset
            if dataset is None:
                num_simulations = len(gdat_files)
                num_observables = len(ssa_result.dtype)
                num_timepoints = len(ssa_result)
                f = h5py.File('%s.hdf5' % filename, 'w')
                dataset = f.create_dataset('data',
                               (num_conditions, num_simulations,
                                num_observables, num_timepoints),
                               dtype='float',
                               chunks=(1, min(SIM_CHUNK_SIZE, num_simulations),
                                       1, min(TIME_CHUNK_SIZE, num_timepoints)),
                               compression=9, shuffle=True)
            # Make sure there's no funny business
            assert len(gdat_files) == num_simulations
            assert len(ssa_result.dtype) == num_observables
            assert len(ssa_result) == num_timepoints
            # Load the data into the appropriate slot in the dataset
            dataset[condition_index, sim_index, :, :] = \
                         ssa_result.view('float').reshape(num_timepoints, -1).T
            print "\rLoaded BNG file %d of %d" % (sim_index+1, num_simulations),
            sys.stdout.flush()
        # (end iteration over all .gdat files in the directory)
        print # Add a newline
    # (end iteration over all data directories)

    # Pickle the dtype with the observables and save in the dataset
    dtype_pck = pickle.dumps(ssa_result.dtype)
    f.create_dataset('dtypes', dtype='uint8', data=map(ord, dtype_pck))
    f.close()
    return dataset

def run_main(jobs):
    """A main method for use in run scripts for job submissions.

    Parases the job index to run from ``sys.argv[1]`` and runs the appropriate
    job from the job list. In addition, creates a subdirectory ``data_%d``,
    where `%d` denotes the job index, and changes the working directory to that
    directory so that SSA simulation results can be output there.

    Parameters
    ----------
    jobs : list of :py:class:`Job` instances
        The list of simulation jobs to run, each one corresponding to a
        different simulation condition.
    """
    # Make sure that we have an argument for the condition index
    if len(sys.argv) <= 1:
        print "The condition index must also be specified."
        sys.exit()
    # Get the job for the condition
    job_index = int(sys.argv[1])
    job = jobs[job_index]
    # Make a directory to put the simulation data in for this condition
    data_dir = os.path.join(os.getcwd(), 'data_%d' % job_index)
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    os.chdir(data_dir)
    # Run the simulation
    job.run_n_cpt(cleanup=False)

if __name__ == '__main__':
    usage_msg =  "Usage:\n"
    usage_msg += "\n"
    usage_msg += "    python simulation.py parse [hdf5_filename] " \
                                                "[files or dirs...]\n"
    usage_msg += "        Parses the BNG simulation results from a list of\n"
    usage_msg += "        directories containing .gdat files and saves the\n"
    usage_msg += "        results as an HDF5 file.\n"
    usage_msg += "\n"
    usage_msg += "    python simulation.py submit [run_script.py] [num_jobs]\n"
    usage_msg += "        For use with LSF. Calls bsub to submit the desired\n"
    usage_msg += "        number of instances of run_script.py. Note that if\n"
    usage_msg += "        run_script.py contains a top-level variable 'jobs'\n"
    usage_msg += "        (which should be a list of instances of a\n"
    usage_msg += "        simulation.Job subclass), the submit script runs\n"
    usage_msg += "        num_jobs jobs for each entry in the job list, for\n"
    usage_msg += "        a total of num_jobs * len(jobs) submitted jobs.\n"

    if len(sys.argv) < 2:
        print usage_msg
        sys.exit()
    # Parse simulation output files
    if sys.argv[1] == 'parse':
        if len(sys.argv) < 4:
            print "Please provide an output filename and a list of files.\n"
            print usage_msg
            sys.exit()
        hdf5_filename = sys.argv[2]
        data_files = sys.argv[3:]
        if hdf5_filename.endswith('.gdat') or os.path.isdir(hdf5_filename):
            print "The argument after 'parse' should be the basename of the " \
                  "HDF5 output file.\n"
            print usage_msg
            sys.exit()
        if len(data_files) == 0:
            print "Please provide a list of files to parse.\n"
            print usage_msg
            sys.exit()
        if os.path.isdir(data_files[0]):
            save_bng_dirs_to_hdf5(data_files, hdf5_filename)
        else:
            print "Please provide a list of .gdat directories."
            print usage_msg
            sys.exit()
    # Submit jobs to LSF
    elif sys.argv[1] == 'submit':
        if len(sys.argv) < 4:
            print usage_msg
            sys.exit()

        run_script = sys.argv[2]
        if not run_script.endswith('.py'):
            print "The run script must be a Python script with .py extension.\n"
            print usage_msg
            sys.exit()

        # Import the run script to see if it has multiple conditions
        import imp
        run_mod = imp.load_source('run_mod', run_script)
        # If the module has a variable 'jobs', we're in a multi-condition
        # situation
        if 'jobs' in run_mod.__dict__.keys():
            num_conditions = len(run_mod.jobs)
        # Otherwise set the num_conditions to 1
        else:
            num_conditions = 1

        # Get the job_name as the output file basename
        if 'job_name' in run_mod.__dict__.keys():
            output_base = run_mod.job_name
        else:
            output_base = run_script.split('.')[0]

        num_jobs = int(sys.argv[3])
        queue = 'short'
        time_limit = '12:00'
        cmd_list = []
        for condition_index in range(num_conditions):
            for job_index in range(num_jobs):
                output_filename = '%s_cond%d_job%d.out' % \
                                  (output_base, condition_index, job_index)

                cmd_list = ['bsub',
                            '-W', time_limit,
                            '-q', queue,
                            '-o', output_filename,
                            'python', run_script,
                            str(condition_index)]
                print ' '.join(cmd_list)
                subprocess.call(cmd_list)
    # Unrecognized command
    else:
        print usage_msg
        sys.exit()


