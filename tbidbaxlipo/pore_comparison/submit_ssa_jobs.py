"""
A script for submitting multiple MCMC jobs to LSF.  Allows fitting multiple
models to multiple combinations of data trajectories with multiple combinations
of model observables.
"""
import subprocess
import numpy as np
import sys

if __name__ == '__main__':
    usage_msg =  "Usage:\n"
    usage_msg += "    python submit_ssa_jobs.py [run_script.py] [num_sims]\n"

    if len(sys.argv) < 3:
        print usage_msg
        sys.exit()
    if not sys.argv[1].endswith('.py'):
        print "The run script must be a Python script with .py extension.\n"
        print usage_msg
        sys.exit()

    run_script = sys.argv[1]
    num_sims = int(sys.argv[2])
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

