Raw data
========

By measurement type
-------------------

.. plot::

    from tbidbaxlipo.plots.nbd_bax_analysis import plot_all
    from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
    plot_all(df, nbd_residues, ['Release', 'NBD'])

By replicate
------------

.. plot::

    from tbidbaxlipo.plots.nbd_bax_analysis import plot_all_by_replicate
    from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
    plot_all_by_replicate(df, nbd_residues, ['Release', 'NBD'])

