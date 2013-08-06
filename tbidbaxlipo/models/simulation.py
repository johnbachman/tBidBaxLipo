from pysb import bng, kappa
import sys
import numpy as np
import pickle
from tbidbaxlipo.models import one_cpt, n_cpt, site_cpt
from pysb.integrate import Solver
import collections

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

    def build(self, module):
        """Virtual method to build the model from the appropriate model Builder.

        Must be implemented by any subclass. Should take a module as object
        (e.g., :py:mod:`tbidbaxlipo.models.one_cpt` or
        :py:mod:`tbidbaxlipo.models.n_cpt`) and return an instance of a Builder
        class after model construction is completed.  """

        raise NotImplementedError()

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
        b = self.build(n_cpt)
        xrecs = []
        for i in range(self.num_sims):
            print "Running BNG simulation %d of %d..." % (i, self.num_sims)
            xrecs.append(bng.run_ssa(b.model, t_end=self.tmax,
                         n_steps=self.n_steps, cleanup=cleanup, output_dir='.'))
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
        b = self.build(site_cpt)
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

        b = self.build(one_cpt)
        # We add one extra point so that the time coordinates of the
        # one_cpt and n_cpt/site_cpt simulations line up
        t = np.linspace(0, self.tmax, self.n_steps + 1)
        s = Solver(b.model, t)
        s.run()
        return (t, s.yobs)

# Operations on stochastic simulation results

def get_observables_values(observable_names, output, timepoint=None):
    """Get the values for a list of observables from a simulation result.

    Useful primarily for getting a list of observables corresponding to a
    given observable distributed across all compartments and then analyzing
    the distribution.

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
    numpy.array
        An array with the simulated values from the output array for the
        list of observables at the given timepoint.
    """

    values = []
    num_timepoints = len(output['time'])
    if timepoint is None:
        timepoint = num_timepoints - 1
    if timepoint < 0 or timepoint >= num_timepoints:
        raise ValueError('Invalid value for the timepoint index.')
    for obs_name in observable_names:
        values.append(output[obs_name][timepoint])
    return np.array(values)

def get_means_across_cpts(observable_names, output, timepoint=None):
    means_list = []
    var_list = []
    for sim_result in output:
        values = get_observables_values(observable_names, sim_result,
                                        timepoint=timepoint)
        means_list.append(np.mean(values))
        var_list.append(np.var(values))
    return (np.array(means_list), np.array(var_list))

def get_frequency_matrix(observable_names, output, timepoint=None):
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
        The first entry is an array of indices giving the values over which the
        frequencies are observed. The second entry is a matrix with shape
        (len(index), len(output)), where (i, j) gives the number of times the
        the value index(i) was observed across the list of observables
        in simulation output(j).
    """

    # The list of counts
    counts_list = []
    for sim_result in output:
        # Get the array with the amounts of each observable in the list
        values = get_observables_values(observable_names, sim_result,
                                        timepoint=timepoint)
        # Convert to a dict of frequencies of each amount
        counts_list.append(collections.Counter(values))

    # Get the set containing the amounts observed across all simulations
    all_keys = set()
    for counts in counts_list:
        all_keys |= set(counts.keys())

    # Convert to a matrix in which the 0th entry represents the frequency
    # of observing the smallest observed amount, and the last entry represents
    # the frequency of observing the largest observed amount
    key_min = min(all_keys)
    key_max = max(all_keys)
    freq_matrix = np.zeros((key_max - key_min + 1, len(counts_list)))
    for i, counts in enumerate(counts_list):
        for key, val in counts.iteritems():
            freq_matrix[key - key_min, i] = val
    # The index runs from the minimum to the maximum amount
    index = np.array(range(int(key_min), int(key_max) + 1))
    return (index, freq_matrix)

def load_bng_files(gdat_files):
    """Load the observable record arrays from BNG simulations.

    Parameters
    ----------
    gdat_files : list of strings
        A list of filenames to be loaded, e.g. as returned by a glob.

    Returns
    -------
    list of numpy.recarrays
        The entries in the list are the results of each simulation, as
        a numpy record array indexed by observable name.
    """

    print "Loading BNG simulation result files..."
    xrecs = []
    # Load all of the simulation results into a list of record arrays
    for gdat_file in gdat_files:
        #print "Loading %s" % gdat_file
        ssa_result = bng._parse_bng_outfile(gdat_file)
        xrecs.append(ssa_result)
    return xrecs

def calculate_mean_and_std(xrecs):
    """Returns the mean and std of the observables from SSA simulations.

    Parameters
    ----------
    xrecs : list of numpy.recarrays
        A list of observable record arrays, for example of the type
        returned by :py:meth:`Job.run_n_cpt`.

    Returns
    -------
    tuple of (numpy.recarray, numpy.recarray)
        The first entry in the tuple is the record array of means
        (across the set of simulations run); the second is the record array
        of standard deviations.
    """

    # Convert the list of record arrays into a matrix
    xall = np.array([x.tolist() for x in xrecs])
    # Create new record arrays with the means and SDs for all of the the
    # observables as the entries
    means = np.recarray(xrecs[0].shape, dtype=xrecs[0].dtype,
                        buf=np.mean(xall, 0))
    stds = np.recarray(xrecs[0].shape, dtype=xrecs[0].dtype,
                        buf=np.std(xall, 0))
    return (means, stds)

def calculate_mean_and_std_from_files(gdat_files):
    """Calculates mean and SD for all observables from .gdat files.

    Parameters
    ----------
    gdat_files : list of strings
        The list of gdat files to load.

    Returns
    -------
    tuple of (numpy.recarray, numpy.recarray)
        The first entry is the record array of means; the second is
        the record array of standard deviations.
    """
    xrecs = load_bng_files(gdat_files)
    return Job.calculate_mean_and_std(xrecs)

def calculate_dye_release_mean_and_std(xrecs, pore_obs_prefix='pores_'):
    # Get the timepoins
    num_timepoints = len(xrecs[0])
    # Get the list of pore observables; makes the assumption that all pore
    # observable names have the same format, with a prefix followed by the
    # compartment identifier, e.g. 'pores_c38'.
    pore_obs_list = [field_name for field_name in xrecs[0].dtype.fields
                     if field_name.startswith(pore_obs_prefix)]
    # Assume that there is a pore observable for every vesicle:
    num_vesicles = len(pore_obs_list)
    # For every simulation result, calculate a dye release vector
    dye_release = np.zeros((len(xrecs), num_timepoints))
    for i, xrec in enumerate(xrecs):
        assert len(xrec) == num_timepoints # all sims must be same length
        # Calculate the fraction of vesicles permeabilized at this timepoint
        for t in range(num_timepoints):
            permeabilized_count = 0
            for pore_obs in pore_obs_list:
                if xrec[pore_obs][t] > 0:
                    permeabilized_count += 1
            dye_release[i, t] = permeabilized_count / float(num_vesicles)
    # Calculate the mean and SD across the matrix
    mean = np.mean(dye_release, axis=0)
    std = np.std(dye_release, axis=0)
    return (mean, std)

if __name__ == '__main__':
    usage_msg =  "Usage:\n"
    usage_msg += "\n"
    usage_msg += "    python simulation.py parse [files]\n"
    usage_msg += "        Parses the BNG simulation results from a list of\n"
    usage_msg += "        .gdat files and saves the results as a list of\n"
    usage_msg += "        record arrays in n_cpt_obs.pck. Also calculates\n"
    usage_msg += "        and saves the mean and SD for all observables in\n"
    usage_msg += "        means.pck and stds.pck, respectively.\n"
    usage_msg += "\n"
    usage_msg += "    python simulation.py submit [run_script.py] [num_jobs]\n"
    usage_msg += "        For use with LSF. Calls bsub to submit the desired\n"
    usage_msg += "        number of instances of run_script.py."

    if len(sys.argv) < 2:
        print usage_msg
        sys.exit()
    # Parse simulation output files
    if sys.argv[1] == 'parse':
        gdat_files = sys.argv[2:]
        if len(gdat_files) == 0:
            print "Please provide a list of files to parse.\n"
            print usage_msg
            sys.exit()
        n_cpt_obs = load_bng_files(gdat_files)
        (means, stds) = calculate_mean_and_std(n_cpt_obs)
        # Pickle the output files
        with open('n_cpt_obs.pck', 'w') as f:
            print "Writing record arrays..."
            pickle.dump(n_cpt_obs, f)
        with open('means.pck', 'w') as f:
            print "Writing means..."
            pickle.dump(means, f)
        with open('stds.pck', 'w') as f:
            print "Writing SDs..."
            pickle.dump(stds, f)
    # Submit jobs to LSF
    elif sys.argv[1] == 'submit':
        if len(sys.argv) < 4:
            print usage_msg
            sys.exit()
        if not sys.argv[1].endswith('.py'):
            print "The run script must be a Python script with .py extension.\n"
            print usage_msg
            sys.exit()

        run_script = sys.argv[2]
        num_sims = int(sys.argv[3])
        queue = 'mini'
        time_limit = '00:10'
        output_base = run_script.split('.')[0]

        cmd_list = []
        for i in range(num_sims):
            output_filename = '%s_%d.out' % (output_base, i)

            cmd_list = ['bsub',
                        '-W', time_limit,
                        '-q', queue,
                        '-o', output_filename,
                        'python', run_script]
            print ' '.join(cmd_list)
            subprocess.call(cmd_list)
    # Unrecognized command
    else:
        print usage_msg
        sys.exit()
