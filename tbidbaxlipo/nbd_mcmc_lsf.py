import sys

num_chains = 20
queue = 'sorger_1d'

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Please specify a fit type (c3, c3c120, or linear)."
    else:
        if sys.argv[1] == 'c3':
            for i in range(0, num_chains):
                cmdstr = "bsub -q %s python nbd_mcmc.py " \
                         "nbd_mcmc_c3_random_initial_values.pck %d nbd_mcmc" % \
                         (queue, i)
                print cmdstr
        elif sys.argv[1] == 'c3c120':
            print "Not implemented yet"
        elif sys.argv[1] == 'linear':
            for i in range(0, num_chains):
                cmdstr = "bsub -q %s python nbd_mcmc.py " \
                         "nbd_mcmc_linear_random_initial_values.pck %d " \
                         "nbd_mcmc_linear" % \
                         (queue, i)
                print cmdstr
        else:
            print "Not a known fit type."
