for i in range(0, 5):
    cmdstr = "bsub -q sorger_1d python nbd_mcmc.py " \
             "nbd_parallel_random_initial_values.pck %d nbd_mcmc" % i
    print cmdstr
