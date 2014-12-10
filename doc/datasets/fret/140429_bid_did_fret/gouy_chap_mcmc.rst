Fits to Gouy-Chapman model
==========================

Model taken from:

Schwarz, G., & Beschiaschvili, G. (1989). Thermodynamic and kinetic studies on
the association of melittin with a phospholipid bilayer. Biochimica et
biophysica acta, 979(1), 82–90.

Fit to data
-----------

First, we load the (previously run) MCMC chain to get our distributions of
parameters:

.. ipython::

    In [1]: from tbidbaxlipo.plots.x140429_Bid_membrane_FRET.gouy_chap_mcmc import *

    In [2]: with open('../results/mcmc/140429_gouy_chap_mcmc.pck') as f:
       ...:     sampler = pickle.load(f)
       ...:

    In [3]: plt.close('all')

    In [3]: plot_chain(sampler.flatchain[0], sampler.lnprobability[0])

We plot the posterior probability of the chains to show that they've
converged, which is indicated by the fact that they are not trending upwards
to better fits:

.. ipython::

    @savefig 140429_gouy_chap_mcmc_1.png
    In [4]: plt.figure('Chains')

Next we plot the model fits to the FRET data. The red curve is the maximum
likelihood fit, while the shaded gray area gives the 95% confidence interval
around the predictions:

.. ipython::

    In [4]: plt.figure('Fits')

    In [5]: plt.savefig('_static/140429_gouy_chap_mcmc_2.png', dpi=300)

.. image:: ../../../_static/140429_gouy_chap_mcmc_2.png
    :width: 400px

Marginal distributions of parameters
------------------------------------

We plot the triangle plot of the parameters:

.. ipython::

    @savefig 140429_gouy_chap_mcmc_3.png width=600
    In [4]: plt.figure(3)

Predictions for other types of binding experiments
--------------------------------------------------

Liposome titration experiment predictions for 20 nM cBid:

.. ipython::

    In [1]: plot_saturation_binding_predictions(sampler.flatchain[0])

    In [2]: plt.figure('Lipobinding')

    In [5]: plt.savefig('_static/140429_gouy_chap_mcmc_5.png', dpi=300)

.. image:: ../../../_static/140429_gouy_chap_mcmc_5.png

Interestingly, the liposome binding experiment looks comparable to the
Bid/liposome binding experiment reported in Shamas-Din et. al, which reported a
KD of ~1 nM liposomes when measured in this fashion.

