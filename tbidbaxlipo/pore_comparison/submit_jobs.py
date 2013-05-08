"""
A script for submitting multiple MCMC jobs to LSF.  Allows fitting multiple
models to multiple combinations of data trajectories with multiple combinations
of model observables.
"""
import subprocess
import numpy as np

# Use this if you want to try fitting every model
# model_names = nbd_mcmc.model_names


if __name__ == '__main__':
    num_sims = 20
    queue = 'short'
    cmd_list = []
    for i in range(0, num_sims):
        output_filename = 'run_n_cpt_job_%d.out' % i

        cmd_list = ['bsub', '-W', '01:00',
                    '-q', queue,
                    '-o', output_filename,
                    'python', '-m', 'tbidbaxlipo.pore_comparison.run_n_cpt']
        print ' '.join(cmd_list)
        subprocess.call(cmd_list)

