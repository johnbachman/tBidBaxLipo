Raw data
========

By measurement type
-------------------

.. plot::

    from tbidbaxlipo.plots.nbd_bax_analysis import plot_all
    from tbidbaxlipo.plots.x141203_Bax_Bax_FRET.preprocess_data \
            import df, nbd_residues
    plot_all(df, nbd_residues, ['Release', 'NBD', 'FRET'], activators=['Bid'],
             replicates=(1,))

By replicate
------------

.. plot::

    from tbidbaxlipo.plots.nbd_bax_analysis import plot_all_by_replicate
    from tbidbaxlipo.plots.x141203_Bax_Bax_FRET.preprocess_data \
            import df, nbd_residues
    plot_all_by_replicate(df, nbd_residues, ['Release', 'NBD', 'FRET'],
                          activators=['Bid'], replicates=(1,))

