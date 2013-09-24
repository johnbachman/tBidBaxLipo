.. _130923_dextran_mDivi1_analogs:

FITC-dextran release with mDivi-1 analogs (9/23/13)
===================================================

The idea behind this experiment was the hypothesis that the partial inhibition
of permeabilization :ref:`observed for the mdivi-1 analogs in liposomes
<130830_mDivi1_analogs>` was due to the fact that the analogs prevented large
pores but permitted pores small enough to release ANTS. In this experiment I
used release of FITC dextrans into a solution containing an
anti-fluorescein antibody (which quenches fluorescein fluorescence upon
binding) as the marker of permeabilization to see if the effects of these
molecules were stronger.

Interestingly, in this assay the drugs showed even **less of an effect** than
in the ANTS assay, though the rank ordering was still the same, with mdivi-1
itself the most potent.

Since the concentration of liposomes here seems likely significantly less than
the concentration I was able to achieve using the ANTS assay, it's possible
that the difference is due to a change in the competitive binding between
mdivi-1 and lipids, and Bax and lipids.

.. plot::

    from tbidbaxlipo.plots.layout_130923 import plot_data
    plot_data()


