bsub -q sorger_1d -o nbc_mcmc_0.out -u bachman@fas.harvard.edu python nbd_linear_mcmc.py nbd_parallel_random_initial_values.pck 0 nbd_mcmc
bsub -q sorger_1d -o nbc_mcmc_1.out -u bachman@fas.harvard.edu python nbd_linear_mcmc.py nbd_parallel_random_initial_values.pck 1 nbd_mcmc
bsub -q sorger_1d -o nbc_mcmc_2.out -u bachman@fas.harvard.edu python nbd_linear_mcmc.py nbd_parallel_random_initial_values.pck 2 nbd_mcmc
bsub -q sorger_1d -o nbc_mcmc_3.out -u bachman@fas.harvard.edu python nbd_linear_mcmc.py nbd_parallel_random_initial_values.pck 3 nbd_mcmc
bsub -q sorger_1d -o nbc_mcmc_4.out -u bachman@fas.harvard.edu python nbd_linear_mcmc.py nbd_parallel_random_initial_values.pck 4 nbd_mcmc
