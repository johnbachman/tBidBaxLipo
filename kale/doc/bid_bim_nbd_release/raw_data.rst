Raw data
========

.. plot::

    from tbidbaxlipo.plots.nbd_bax_analysis import plot_all
    from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
    plot_all(df, nbd_residues, ['Release', 'NBD'])
