for i in range(0, 5):
    cmdstr = "bsub -q sorger_1d -o nbc_mcmc_%d.out python nbd_linear_mcmc.py " \
             "nbd_parallel_random_initial_values.pck %d nbd_mcmc" % (i, i)
    print cmdstr
