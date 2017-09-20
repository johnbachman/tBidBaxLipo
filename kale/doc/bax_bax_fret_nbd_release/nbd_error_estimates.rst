Estimates of error in NBD curves
================================

.. plot::

    from tbidbaxlipo.data.parse_bax_bax_fret_nbd_release import df, nbd_residues
    from tbidbaxlipo.plots.nbd_bax_analysis import plot_nbd_error_estimates
    plot_nbd_error_estimates(df, nbd_residues, last_n_pts=80, fit_type='cubic')

