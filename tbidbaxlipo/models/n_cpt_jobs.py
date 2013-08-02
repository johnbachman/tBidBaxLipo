from pysb import bng
import sys
import numpy as np
import pickle

class Job(object):
    """Class for storing aspects of a SSA simulation job.

    Parameters
    ----------
    model : instance of pysb.core.Model
        The model to simulate.
    tmax : number
        The duration of the simulation.
    n_steps : int
        The number of simulation data points to save and return.
    n_sims : int
        The number of simulations to run as part of this job.
    """
    def __init__(self, model, tmax, n_steps, n_sims):
        self.model = model
        self.tmax = tmax
        self.n_steps = n_steps
        self.n_sims = n_sims

    def run(self):
        for i in range(self.n_sims):
            bng.run_ssa(self.model, t_end=self.tmax, n_steps=self.n_steps,
                        cleanup=False, output_dir='.')

def calculate_mean_and_std_from_files(gdat_files):
    """Calculates mean and SD for all observables from .gdat files.

    The results are saved as pickled output files.

    Parameters
    ----------
    gdat_files : list of strings
        The list of gdat files to load.

    Returns
    -------
    list of numpy record arrays
        The first entry is the record array of means; the second is
        the record array of standard deviations.
    """

    print "Calculating mean and SD for observables from .gdat files..."
    xrecs = []
    dr_all = []
    # Load all of the simulation results into a list of record arrays
    for gdat_file in gdat_files:
        print "Loading %s" % gdat_file
        ssa_result = bng._parse_bng_outfile(gdat_file)
        xrecs.append(ssa_result)
        #dr_all.append(get_dye_release(model, 'pores', ssa_result))
    # Convert the list of record arrays into a matrix
    xall = np.array([x.tolist() for x in xrecs])
    # Create new record arrays with the means and SDs for all of the the
    # observables as the entries
    means = np.recarray(xrecs[0].shape, dtype=xrecs[0].dtype,
                        buf=np.mean(xall, 0))
    stds = np.recarray(xrecs[0].shape, dtype=xrecs[0].dtype,
                        buf=np.std(xall, 0))
    return [means, stds]
    #with open('dr.pck', 'w') as f:
    #    pickle.dump(dr_all, f)

if __name__ == '__main__':
    usage_msg =  "Usage:\n"
    usage_msg += "\n"
    usage_msg += "    python ssa_job.py parse [files]\n"
    usage_msg += "        Calculates the mean and SD for all observables\n"
    usage_msg += "        from a list of .gdat files and saves the results in\n"
    usage_msg += "        two pickle files, one for the mean and one for SD.\n"
    usage_msg += "\n"
    usage_msg += "    python ssa_job.py submit [run_script.py] [num_jobs]\n"
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
        (means, stds) = calculate_mean_and_std_from_files(gdat_files)
        # Pickle the output files
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
