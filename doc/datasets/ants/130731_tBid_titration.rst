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
permeabilization--and this even though the concentration of liposomes was 
24 nM. Does this perhaps suggest that tBid can indeed unbind from liposomes to
go on to permeabilize other liposomes? But then why was there no unbinding signal
in the octet experiment? Perhaps tBid can only unbind when Bax is present?

.. plot::

    from tbidbaxlipo.plots.layout_130731 import main
    main()


