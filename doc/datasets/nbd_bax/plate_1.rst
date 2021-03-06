NBD-Bax Kinetics, Plate Reader
==============================

A dataset from Justin Kale in which multiple mutants were assayed at the same
time in a plate reader. The dataset includes 15 mutants, with four replicates
for each. Replicates 0 and 1 were performed together on the same day;
replicates 2 and 3 were performed together on a different day.

Unlike the fluorimeter dataset, in these experiments the cBid was first
equilibrated with the liposomes and then the Bax was added afterwards. The
values given here are the raw fluorescence values (not normalized to
fold-change over background).

Raw data
--------

.. plot::

    from tbidbaxlipo.data.nbd_plate_data import data
    from matplotlib import pyplot as plt

    def plot_raw(data, display=False):
        # For each mutant, make a plot of all replicates
        nbd_names = data.columns.levels[0]
        replicates = data.columns.levels[1]
        for nbd_name in nbd_names:
            plt.figure()
            for replicate in replicates:
                tc = data[(nbd_name, replicate)]
                plt.plot(tc[:, 'TIME'], tc[:, 'VALUE'],
                         label='%s rep. %d' % (nbd_name, replicate))
                plt.legend(loc='lower right')
                plt.xlabel('Time (sec)')
                plt.ylabel('Fluorescence (AU)')
                plt.title(nbd_name)

    plot_raw(data)

Notes on Unimolecular Multi-Conformer Fits, 7/16/13
---------------------------------------------------

Notes on fits to the data performed by MCMC using the "unimolecular
multi-conformer" model. An index of the fits, performed on 7/16/13, can be
found on the `Sorger lab web site
<http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/index.html>`_.

c5_
    Has a pretty good fit just with two conformations, though the error is
    reduced most sharply in the switch from 2 to 3 conformations. However, the
    improvement in the autocorrelation plot (and in the Durbin-Watson
    statistic) between the conformations is not especially dramatic.
    **Conclusion:**  Not definitive; all have fairly good fits.

c15_
    Even though the fit with two conformations looks decent, there is a very
    sharp reduction in error between 2 and 3, and an even sharper improvement
    in the residuals in going from 2 to 3. However, there is only marginal
    improvement in either in going above 3 conformations.  **Conclusion:**
    Strong evidence for 3 conformations.

c36_
    Fits are not great with any number of conformations; all of the ACF plots
    show systematic error. None of the models can fit the big bump in
    replicates 0 and 1. For some of the replicates the 4 conformation model is
    able to fit the initial transient. Data is noisy
    **Conclusion:** Not definitive, all have fairly poor fits.

c40_
    Data is very noisy; all fits are fairly poor.  **Conclusion:** Not
    definitive.

c47_
    Fits at 2 conformations look good, but there is a clear indication from the
    ACF plots that switching to 3 or more conformations eliminates systematic
    error. There is minimal improvement from additional conformations.
    **Conclusion:** Strong evidence for 3 conformations.

c54_
    Interestingly, the c54 mutant is one of the few that indicates 4
    conformations, at least for the replicates 1, 2, and 3 (not 0). The data
    has a sharp initial transient followed by a drop and a slow rise. There are
    sharp drops in error and sharp improvements in the ACF plot between both 2
    and 3, and 3 and 4 conformations.  **Conclusion:** Fairly strong evidence
    for 4 conformations.

c62_
    For replicates 0 and 1, the switch from 2 to 3 conformations seems to
    dramatically improve both the fit and the ACF. However, the improvement is
    more moderate for replicates 2 and 3. The data is fairly noisy and could be
    improved by fluorimeter-based measurements.  **Conclusions:** Moderate
    evidence for 3 conformations.

c68_
    Marked improvement from 2 to 3, and a subtler but still seemingly
    significant improvement in fit from 3 to 4.  **Conclusions:** Strong
    evidence for at least 3 states; weaker evidence for 4 conformations.

c79_
    For replicate 0, the big bump around 2000 seconds makes the error for all
    fits seem systematic. However, the biggest improvement in fit is between 2
    and 3 conformations. For replicates 1 and 2, the fits seem to definitively
    suggest 3 conformations and no more. For replicate 3, at least 3
    conformations are required, with a fairly slight improvement with 4.
    **Conclusion:** Strong evidence for 3 conformations; weak evidence for 4.

c120_
    Fit at 2 conformations is already very good; adding additional
    conformations marginally improves it each time. This is fairly consistent
    with the previous fluorimeter-based data, in which c120 was the mutant that
    had quite a good fit with a single exponential (though the fit still
    improved with a double exponential). **Conclusion:** No strong evidence for
    more than 2 conformations.

c122_
    **Conclusion:** Strong evidence for 3 conformations, and no more.

c126_
    Much as with c120_, the fit at 2 conformations is already fairly good, with
    marginal improvements with additional conformations. The shoulder in the
    trace for replicate 0 (and to some extent replicate 1 also) means that the
    error appears systematic for all of these. For replicates 2 and 3, there is
    somewhat of an improvement from 2 to 3, but this might be mitigated if the
    fluorescence of the initial state were allowed to deviate slightly from the
    initial fluorescence value. **Conclusion:** No strong evidence for more
    than 2 conformations.

c175_
    Much as with c120_ and c126_, the fit at 2 conformations is fairly good but
    improves with additional conformations. While the error is seemingly
    systematic for all conformations, by far the most substantial decrease in
    error comes from the change from 2 to 3 conformations. This seems to be due
    at least in part to the need to fit a sharp initial transient occuring in
    the first few timepoints.  **Conclusion:** Weak evidence for 3
    conformations (and not more).

c179_
    A biphasic curve with a sharp drop and then a slow rise. Absolutely cannot
    be fit by a 2 conformation model. Fits with 3 or more are quite good,
    though the improvement with more than 3 is fairly minimal. **Conclusion:**
    Very strong evidence for 3 conformations.

c188_
    Though the curve appears fairly like a straightforward exponential, there
    is a very definitive improvement in fit between 2 and three conformations
    indicated both by the error and the ACF plots. Moreover, there is very
    minimal improvement with more than three conformations. **Conclusion:**
    Very strong evidence for 3 conformations.

.. _c5: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c5/index.html
.. _c15: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c15/index.html
.. _c36: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c36/index.html
.. _c40: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c40/index.html
.. _c47: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c47/index.html
.. _c54: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c54/index.html
.. _c62: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c62/index.html
.. _c68: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c68/index.html
.. _c79: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c79/index.html
.. _c120: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c120/index.html
.. _c122: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c122/index.html
.. _c126: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c126/index.html
.. _c175: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c175/index.html
.. _c179: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c179/index.html
.. _c188: http://sorger.med.harvard.edu/data/bachman/130716_nbd_plate_fits/c188/index.html

Summary
^^^^^^^

All conclusions based solely on this dataset.

Mutants without strong evidence for more than 2 conformations: c5_,  c120_, c126_

Evidence for 3 conformations: c15_, c47_, c62_, (c68_, c79_), c122_, c175_, c179_, c188_

Evidence for 4 conformations: c54_, c68_, c79_

Data noisy, inconclusive: c36_, c40_

Locations of the Different Cysteines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As indicated by the Bax solution structure, PDBID 1F16:

* c3, c5: unstructured N-terminus
* c15: N-term end of a1
* c36: C-term end of a1
* c40, c47: unstructured region between a1 and a2
* c54, (c62, c68): a2 (BH3)
* c79: a3
* c120, c122, c126: a5
* c175, c179, c188: a9

Notes on Unimolecular Multi-Conformer Fits by PT, 7/19/13
---------------------------------------------------------

In this case I ran much longer fits (500000 steps) by parallel tempering to
ensure convergence, with the goal of being able to compare the parameter values
between them. For all residues only one replicate was fit; obviously I'll need
to fit the replicates to know the within-mutant vs. between mutant variability
in the fitted kinetic constants.

c5

    Minimal value of the intermediate state. The time constant for the
    C0 to C1 transition at ML was 5.03e-4 (with log10 of -3.29).

c15

    Bimodal, with C0 to C1 at -3.2 (estimated from marginals plot by eye)
    or -2.6 (by eye). Bigger peak at -3.2. For C1 to C2, big peak at ~-2.6,
    narrower but also big peak at -3.2.

c175

    ML fits has C0-C1 at 6.42e-4 (-3.2) (C1 scaling 25177, seemingly an
    unrealistically high value), C1-C2 at 1.97e-2 (C2 scaling 11064). Marginals
    are bimodal, but the dominant mode of C0-C1 is around -3.2. I wonder if
    it's possible that a different change in prior knowledge would suggest that
    the initial transition is what's faster? Perhaps the initial bump in the
    curve is due to the displacement of alpha 9 by the attacking BH3, for
    example.

c179

    Biphasic curve. There is a slight bump in the first few points that makes
    me wonder if there's a signature from the initial displacement. However, in
    general C0-C1 rate has a definitive peak at around -3.5 (0.00037), and
    C1-C2 has a high peak around -2.8. It is interesting that the initial rates
    are in many cases the slower ones; consistent with notion that Bax
    insertion is the rate limiting step? Also interesting is that the C1
    scaling value is totally non-identifiable--there is a lot of leeway in the
    precise value of the C1 fluorescence to get the initial downward transient
    to be fit. ML values: C0-C1 0.00037, C1-C2 0.0017.

c188

    Bimodal fits; interestingly there appear to be two modes, one in which the
    C1 scaling is less then the C2 scaling, one in which it is more than the C2
    scaling. The more robust one (the one that appears at higher temperatures)
    is the one where C0-C1 is around -3.1, C1 scaling is 3.6 (lower than
    final value of ~3.8), C1-C2 similar (3.1).
