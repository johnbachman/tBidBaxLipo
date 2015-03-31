Fits to two-conformation model
================================

In the bar plots of the fitted k1 values, the error bars are the
standard deviations of the best fit values determined from the fit,
while the three bars are fits for replicates 1, 2, and 3.

Bid
---

.. plot::

    from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
    from tbidbaxlipo.plots.nbd_bax_analysis import plot_2conf_fits
    plot_2conf_fits(df, nbd_residues, 'Bid')

Bim
---

.. plot::

    from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
    from tbidbaxlipo.plots.nbd_bax_analysis import plot_2conf_fits
    plot_2conf_fits(df, nbd_residues, 'Bim')


