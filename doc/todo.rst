To-Do list
==========

Near-term todos
---------------

* Fix package structure
* Have model construction insert a stringable name into the builder obj

* Run fits at different temperatures
* Add a plot for how posterior changes for chains, and add legend to link to
  parameter traces as well

* Sphinx report for adding plots to docs?
* Would be cool if make the docs also ran the various fitting/plotting
  functions in nbd_analysis and grid_analysis, and stuck them in .html files
  and added them as a section of the docs

* Add an overview page to the whole repo (intro.rst, and maybe a README.rst)
* Add an overview page in /models/index.rst
* Very brief overview page for analysis.rst

* Write nbd_mcmc to work with linear/parallel as well as mechanistic models
* Tests to make sure that nbd_mcmc will run

* Make colors of data and timecourses match in fit_plotting_function
* Sphinx plugin to show .evidence field of a function as a part of the
  documentation?
* Problem with scaling factors and normalization by Bax_0 or tBid_0 when amount
  of tBid0/Bax0 changes in an "experiment"
* Formulate concrete plan for how to integrate residue sequence into
  mechanistic model--perhaps do simultaneous MOMP assay?
* Understand permeabilization in bulk model vs. in stochastic model
* Science priority--design an experiment to distinguish models for c62
* Use data that has not been normalized
* Need topology test for what the nbd_observable(s) are
* Need ability to run sim, measure Bax signal at membranes, and then add tBid,
  and measure signal again.
* Better priors for rates?

Convergence issues
------------------
* Convergence of sample problem:
    - Hessian blow up? Is it because of np.inf?
    - Better if no annealing phase? Just calculate hessian at time 0 and go

* MCMC on linear model
    - Do fit of multiple chains on the simple linear model, see how well it
      converges

* Write function to restart chains:
    - need ability to load a chain
    - reattach likelihood function to chain (or perhaps better, put likelihood
      in another file)
    - start running again (same random seed?)
    - Write function to pool chains and summarize distributions
    - Need ability to batch restart a group of chains to re-anneal and repeat,
      or to run using a certain interval as their covariance matrix

* PyMC:
    - Figure out why PyMC is getting stuck (may require visualization). Step
      size??  (Tried setting proposal SD for scaling to be low; tried
      AdaptiveMetropolis, neither worked).
    - Figure out why/whether PyMC is slow
    - Look up literature on other adaptive approaches?

Mechanistic questions to address
--------------------------------

* Can we establish whether tBid activates Bax via a rear site or by a
  BH3:groove interaction from kinetic studies?
* Can we determine whether the signal from c62 is from a protein-protein or a
  protein-lipid interface?

* Problem -- handling of 2D diffusion?
* Look at effect of additional liposomes as an experiment (as proteins/liposome
  ratio goes down, how does this affect result?)

* Look at large polymer scenario

Features needed
---------------

* Need ability to match a ComplexPattern against another ComplexPattern to
  do topology reporting
* Thin chains as they go
* When building MCMC sets, need to make sure that the model is the same for all
* Rewrite nose tests for convergence_criterion (or maybe just do away with
  convergence_criterion tests and keep others)
* Ability extract a covariance matrix from a chain or set of chains, or a set
  of positions; ability to load a covariance matrix into mcmc.options for a
  static run
* observables that are functions of other observables, with expressions that
  can include parameters--this would solve the problem of having to separately
  normalize model output to plot alongside data

To code review with Jeremy
--------------------------

* fit_plotting_function (should move fit_plotting_function to bayessb.MCMC?)
* __init__.py : prune, get_mixed_accepts
* multichain.MCMCSet; way pruning and pooling is handled. Possible to create a
  single pooled MCMC object instead of a list of positions?
* convergence takes MCMCSet object now
* report.convergence_criterion wraps convergence.convergence_criterion
* tbidbaxlipo.nbd_mcmc_pysb.import_chain_groups
* Report takes dict of lists rather than list of MCMCSets (so that all mcmcsets
  don't have to be in memory at the same time)
* Use of multiprocessing/MPI for report generation?
* All reporters now take MCMCSet object (though this can contain one chain)
* Way of passing num_samples to reporters from report? Or set in
  module-by-module basis?
* Inconsistent nomenclature (chain/mcmc) in tbid, bayessb, pysb.report
* pattern in which mcmc objects are subclassed to generate specialized
  subclasses which can be pickled to rebuild state of likelihood_fn, step_fn,
  etc.

Other todos
-----------

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

