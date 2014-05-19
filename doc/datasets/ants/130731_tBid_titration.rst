tBid Titration (7/31/13)
========================

The idea behind this experiment was to see if excessive amounts of tBid could
actually slow down Bax pore formation activity. I used two concentrations of
Bax, 63 and 252 nM, and a 10-fold dilution series for tBid, from 2.3 uM down to
0.24 nM.

.. note:: Mistake in liposome concentration!

    I'm fairly sure that in preparing this experiment I accidentally pipetted
    in 50 uL of liposomes into the final reactions instead of the usual 10 uL.
    This is because I didn't dilute the liposome solution in buffer before
    mixing it with the tBid--I added the tBid to the undiluted liposome
    fraction. The maximal fluorescence values are consistent with being about 5
    times the usual max I see for 10 uL lipid (which is around 70,000 counts).
    So this means that the liposome concentration is somewhere around 25 nM.

First and foremost, there is no evidence in this experiment that high tBid can
inhibit Bax. On the contrary, one sees increases in permeabilization kinetics
even when going from 235 nM to 2.35 uM. It begs the question of the
structure/function of the apparent tBid/Bax binding.

Second, it is interesting that 0.24 nM has virtually no activity, whereas the
2.4 nM tBid has a strong effect, reaching about 50% of the maximal level of
permeabilization--and this even though the concentration of liposomes was 24
nM. Does this perhaps suggest that tBid can indeed unbind from liposomes to go
on to permeabilize other liposomes? But then why was there no unbinding signal
in the octet experiment? Perhaps tBid can only unbind when Bax is present?

.. plot::
    :context:

    from tbidbaxlipo.plots.layout_130731 import *
    plot_data()

Two-exponential fits
--------------------

For 63 nM Bax only:

.. plot::
    :context:

    plt.close('all')
    bax63_layout = extract(bax63_wells, layout)
    (fmax_arr, k1_arr, k2_arr, conc_list) = \
          titration_fits.plot_two_exp_fits(norm_wells, bax63_layout,
                                           num_reps=1, conc_str_index=4,
                                           plot=True)
    plt.figure()
    plt.semilogx(conc_list, fmax_arr[0], marker='o')
    plt.xlabel('[cBid] (nM)')
    plt.ylabel('$F_{max}$')

    plt.figure()
    plt.semilogx(conc_list, k1_arr[0], marker='o')
    plt.xlabel('[cBid] (nM)')
    plt.ylabel('$k_1$')

    plt.figure()
    plt.semilogx(conc_list, k2_arr[0], marker='o')
    plt.xlabel('[cBid] (nM)')
    plt.ylabel('$k_2$')

For 252 nM Bax only:

.. plot::
    :context:

    plt.close('all')
    bax252_layout = extract(bax252_wells, layout)
    (fmax_arr, k1_arr, k2_arr, conc_list) = \
          titration_fits.plot_two_exp_fits(norm_wells, bax252_layout,
                                           num_reps=1, conc_str_index=4,
                                           plot=True)
    plt.figure()
    plt.semilogx(conc_list, fmax_arr[0], marker='o')
    plt.xlabel('[cBid] (nM)')
    plt.ylabel('$F_{max}$')

    plt.figure()
    plt.semilogx(conc_list, k1_arr[0], marker='o')
    plt.xlabel('[cBid] (nM)')
    plt.ylabel('$k_1$')

    plt.figure()
    plt.semilogx(conc_list, k2_arr[0], marker='o')
    plt.xlabel('[cBid] (nM)')
    plt.ylabel('$k_2$')

Discussion
----------

The conclusion appears to be that the Fmax value is set by the amount of Bax,
and is relatively insensitive (though not completely insensitive, as seen in
the cBid titration with 252 nM Bax) to the amount of cBid.

In addition, the kinetic increase (in k1) appears to scale roughly log-linearly
with the amount of Bid, with the increase in rate proportional to the
fold-change in the amount of cBid.

One thing that is interesting is that there is no evidence of saturation of the
cBid activity--that is, if the addition of cBid above some stoichiometric level
leads to a plateau in the amount of cBid at liposomes, one might expect that
the Bax activation kinetics might also plateau. However, this appears to not be
the case, with increases in Bax permeabilization activity that increases even
with Bid concentrations of 2 uM!

However, it's important to note that these fits do not take into account the
basal activity of the Bax with no cBid added. If this were taken into account,
the "threshold" level of cBid required to get appreciable Bax activity (above
background) would be more apparent.

