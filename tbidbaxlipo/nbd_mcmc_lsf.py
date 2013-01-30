import sys
import numpy as np
from pysb.integrate import odesolve
import subprocess
import time

# The calls to odesolve are done to insure that the chains don't crash on each other
# trying to run scipy.weave at the same time.

num_chains = 20
queue = 'sorger_1d'

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Please specify a fit type (c3, c3c120, parallel or linear)."
    else:
        # Fit parallel model with c3 only
        if sys.argv[1] == 'c3':
            for i in range(0, num_chains):
                cmd_list = ["bsub", "-q", queue, "python", "nbd_mcmc_c3.py",
                            "nbd_mcmc_c3_random_initial_values.pck",
                            str(i), "nbd_mcmc_c3"]
                print ' '.join(cmd_list)
                subprocess.call(cmd_list)

        # Fit parallel model with c3 and c120 only
        elif sys.argv[1] == 'c3c120':
            for i in range(0, num_chains):
                cmd_list = ["bsub", "-q", queue, "python", "nbd_mcmc_c3c120.py",
                            "nbd_mcmc_c3c120_random_initial_values.pck",
                            str(i), "nbd_mcmc_c3c120"]
                print ' '.join(cmd_list)
                subprocess.call(cmd_list)

        # Fit parallel model
        elif sys.argv[1] == 'parallel':
            for i in range(0, num_chains):
                cmd_list = ["bsub", "-q", queue, "python", "nbd_mcmc.py",
                            "nbd_mcmc_parallel_random_initial_values.pck",
                            str(i), "nbd_mcmc_parallel"]
                print ' '.join(cmd_list)
                subprocess.call(cmd_list)

        # Fit linear model
        elif sys.argv[1] == 'linear':
            for i in range(0, num_chains):
                cmd_list = ["bsub", "-q", queue, "python", "nbd_mcmc.py",
                            "nbd_mcmc_linear_random_initial_values.pck",
                            str(i), "nbd_mcmc_linear"]
                print ' '.join(cmd_list)
                subprocess.call(cmd_list)

        # Error
        else:
            print "Not a known fit type."
