import sys
import numpy as np
from pysb.integrate import odesolve

# The calls to odesolve are done to insure that the chains don't crash on each other
# trying to run scipy.weave at the same time.

num_chains = 20
queue = 'sorger_1d'
t = np.linspace(0, 60*60*3, 100)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Please specify a fit type (c3, c3c120, or linear)."
    else:
        if sys.argv[1] == 'c3':
            from nbd_parallel_model import model
            x = odesolve(model, t)            
            for i in range(0, num_chains):
                cmdstr = "bsub -q %s python nbd_mcmc.py " \
                         "nbd_mcmc_c3_random_initial_values.pck %d nbd_mcmc" % \
                         (queue, i)
                print cmdstr
        elif sys.argv[1] == 'c3c120':
            from nbd_parallel_model import model
            x = odesolve(model, t)            
            for i in range(0, num_chains):
                cmdstr = "bsub -q %s python nbd_mcmc_c3c120.py " \
                         "nbd_mcmc_c3c120_random_initial_values.pck %d " \
                         "nbd_mcmc_c3c120" % \
                         (queue, i)
                print cmdstr
        elif sys.argv[1] == 'linear':
            from nbd_linear_model import model
            x = odesolve(model, t)            
            for i in range(0, num_chains):
                cmdstr = "bsub -q %s python nbd_mcmc.py " \
                         "nbd_linear_random_initial_values.pck %d " \
                         "nbd_mcmc_linear" % \
                         (queue, i)
                print cmdstr
        else:
            print "Not a known fit type."
