Data after outlier removal
==========================

By measurement type
-------------------

.. plot::

    from tbidbaxlipo.plots.nbd_bax_analysis import plot_all
    from tbidbaxlipo.plots.bax_bax_fret_nbd_release.preprocess_data \
            import df_pre, nbd_residues
    plot_all(df_pre, nbd_residues, ['Release', 'NBD', 'FRET'])

By replicate
------------

.. plot::

    from tbidbaxlipo.plots.nbd_bax_analysis import plot_all_by_replicate
    from tbidbaxlipo.plots.bax_bax_fret_nbd_release.preprocess_data \
            import df_pre, nbd_residues
    plot_all_by_replicate(df_pre, nbd_residues, ['Release', 'NBD', 'FRET'])

