Fits to competitive binding equation
====================================

Fitting equation taken from:

Wang, Zhi-Xin, "An exact mathematical expression for describing competitive
binding of two different ligands to a protein molecule." FEBS Letters 360
(1995) 111-114.

.. ipython::

    In [1]: from tbidbaxlipo.plots.x140429_Bid_membrane_FRET.exact_comp_bind_mcmc import *

    In [2]: sampler = run_mcmc()

    In [3]: plot_chain(sampler)

    @savefig exact_comp_bind_mcmc_1.png
    In [4]: plt.figure('Chains')

    @savefig exact_comp_bind_mcmc_2.png
    In [4]: plt.figure('Fits')

