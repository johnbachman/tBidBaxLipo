Estimates of error
==================

NBD curves
----------

.. plot::

    from tbidbaxlipo.plots.x141203_Bax_Bax_FRET.preprocess_data \
            import df, nbd_residues
    from tbidbaxlipo.plots.nbd_bax_analysis import plot_nbd_error_estimates
    plot_nbd_error_estimates(df, nbd_residues, dtype='NBD', activators=['Bid'],
                             replicates=(1,), last_n_pts=80, fit_type='cubic')

After removing outliers
~~~~~~~~~~~~~~~~~~~~~~~

.. plot::

    from tbidbaxlipo.plots.x141203_Bax_Bax_FRET.preprocess_data \
            import df_pre, nbd_residues
    from tbidbaxlipo.plots.nbd_bax_analysis import plot_nbd_error_estimates
    plot_nbd_error_estimates(df_pre, nbd_residues, dtype='NBD',
                             activators=['Bid'], replicates=(1,),
                             last_n_pts=80, fit_type='cubic')


FRET curves
-----------

.. plot::

    from tbidbaxlipo.plots.x141203_Bax_Bax_FRET.preprocess_data \
            import df, nbd_residues
    from tbidbaxlipo.plots.nbd_bax_analysis import plot_nbd_error_estimates
    plot_nbd_error_estimates(df, nbd_residues, dtype='FRET',
                             activators=['Bid'], replicates=(1,),
                             last_n_pts=80, fit_type='cubic')

After removing outliers
~~~~~~~~~~~~~~~~~~~~~~~

.. plot::

    from tbidbaxlipo.plots.x141203_Bax_Bax_FRET.preprocess_data \
            import df_pre, nbd_residues
    from tbidbaxlipo.plots.nbd_bax_analysis import plot_nbd_error_estimates
    plot_nbd_error_estimates(df_pre, nbd_residues, dtype='FRET',
                             activators=['Bid'], replicates=(1,),
                             last_n_pts=80, fit_type='cubic')



