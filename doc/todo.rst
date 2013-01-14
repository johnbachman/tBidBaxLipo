ToDo List for Development
=========================

.. todo:: Fix .rst file for nbd_mcmc (now nbd_linear_mcmc)

.. todo:: Plot the nbd_parallel_model against the normalize_fit data to make
sure that the parameters derived from the single_exp fit in nbd_analysis.py
produce the same results!

.. todo:: Come up with a strategy for how to sample from starting distributions for NBD MCMC

   One approach would be generate a list of random numbers, write them to a
   file (text or pickled), and then have the "dispatcher" script pass an index
   number to the code to be run, which would then retrieve the appropriate
   input parameters from the matrix.

.. todolist::
