Fits to competitive binding equation
====================================

Fitting equation taken from:

Wang, Zhi-Xin, "An exact mathematical expression for describing competitive
binding of two different ligands to a protein molecule." FEBS Letters 360
(1995) 111-114.

Fit to data
-----------

First, we load the (previously run) MCMC chain to get our distributions of
parameters:

.. ipython::

    In [1]: from tbidbaxlipo.plots.x140429_Bid_membrane_FRET.exact_comp_bind_mcmc import *

    In [1]: from os.path import join

    In [1]: import tbidbaxlipo.plots.x140429_Bid_membrane_FRET as bmf

    In [2]: with open(join(bmf.__path__[0], '140429_exact_comp_bind.mcmc')) as f:
       ...:     sampler = pickle.load(f)
       ...:

    In [3]: plt.close('all')

    In [3]: plot_chain(sampler.flatchain[0], sampler.lnprobability[0])

We plot the posterior probability of the chains to show that they've
converged, which is indicated by the fact that they are not trending upwards
to better fits:

.. ipython::

    @savefig 140429_exact_comp_bind_mcmc_1.png
    In [4]: plt.figure('Chains')

Next we plot the model fits to the FRET data. The red curve is the maximum
likelihood fit, while the shaded gray area gives the 95% confidence interval
around the predictions:

.. ipython::

    In [4]: plt.figure('Fits')

    In [5]: plt.savefig('_static/140429_exact_comp_bind_mcmc_2.png', dpi=300)

.. image:: ../../../_static/140429_exact_comp_bind_mcmc_2.png
    :width: 400px

Marginal distributions of parameters
------------------------------------

We plot the triangle plot of the parameters:

.. ipython::

    @savefig 140429_exact_comp_bind_mcmc_3.png width=600
    In [4]: plt.figure(3)

A few observations about these marginals. The marginal distribution of the KD
for binding between cBid and the hypothetical lipid binding site is not constrained with respect to its lower bound, but has a well-defined upper bound, as
shown by the histogram. The maximum value observed by MCMC is:

.. ipython::

    # Max value for Bid/lipo-site KD in nanomolar
    In [3]: print 10 ** np.max(sampler.flatchain[0,:,0])

The notion of a percentile statistic for this upper bound (e.g., a 95%
confidence interval) is not especially meaningful since the 95% interval will
be dependent on the lower bound used for the prior, in this case 0.1 picomolar.
If we had an even smaller lower bound for the KD, this would fill up more area
in the marginal distribution, and bring down the estimate of the 95% upper
bound.

The other relevant parameter, the number of binding sites in the experiment, is
constrained to a defined range. We can calculate the mean and 95%
confidence interval:

.. ipython::

    # Mean [Lipo-sites], in nanomolar
    In [18]: mean_lipo_site_conc = 10 ** np.mean(sampler.flatchain[0,:,1])

    In [19]: print mean_lipo_site_conc

    #  95% conf interval for [Lipo-sites], in nanomolar
    In [18]: lipo_site_concs_95 = 10 ** np.percentile(sampler.flatchain[0,:,1], [2.5, 97.5])

    In [19]: print lipo_site_concs_95

Given the initial concentration of liposomes of 0.1 mg/mL ~= 1.55 nM, this gives us a confidence interval on the number of binding sites per liposome:

.. ipython::

    # Mean sites per liposome
    In [19]: print mean_lipo_site_conc / 1.55

    # Sites per liposome, 95% confidence interval
    In [19]: print lipo_site_concs_95 / 1.55

The resulting values are strikingly low, suggesting a maximum of ~6 Bids
per liposome.

Predictions for other types of binding experiments
--------------------------------------------------

Finally, we use the parameter distributions from the to fit to make predictions
about two different types of experiment: a saturation binding experiment in
which a fluorescent donor is titrated (rather than unlabeled competitor) and a
liposome titration.

First, the Bid titration experiment predictions, for 0.1 mg/mL ~= 1.55 nM liposomes. Here the red line is the mean prediction, and the shaded area is the 95% confidence interval:

.. ipython::

    In [1]: plot_saturation_binding_predictions(sampler.flatchain[0])

    In [2]: plt.figure('Satbinding')

    In [5]: plt.savefig('_static/140429_exact_comp_bind_mcmc_4.png', dpi=300)

.. image:: ../../../_static/140429_exact_comp_bind_mcmc_4.png

This plot highlights the uncertainty about the fraction of cBid bound at low
concentrations, where liposome binding sites are not limiting. To get a better
estimate of this value, it should be possible to increase the concentration of
liposomes (perhaps by 5-fold) and perform the binding experiment so that the
plateau in binding for low cBid concentrations could be observed.

Next, the liposome titration experiment predictions for 20 nM cBid:

.. ipython::

    In [2]: plt.figure('Lipobinding')

    In [5]: plt.savefig('_static/140429_exact_comp_bind_mcmc_5.png', dpi=300)

.. image:: ../../../_static/140429_exact_comp_bind_mcmc_5.png

Interestingly, the liposome binding experiment looks comparable to the
Bid/liposome binding experiment reported in Shamas-Din et. al, which reported a
KD of ~1 nM liposomes when measured in this fashion. In this case the apparent
KD would appear to be a bit higher, perhaps ~3.5 nM on average, but with 1 nM
close to the 95% confidence interval.

