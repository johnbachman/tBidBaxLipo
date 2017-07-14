Raw data
========

By measurement type
-------------------

.. plot::

    from tbidbaxlipo.plots.nbd_bax_analysis import plot_all
    from tbidbaxlipo.data.parse_bax_bax_fret_nbd_release import df, nbd_residues
    plot_all(df, nbd_residues, ['Release', 'NBD', 'FRET'], activators=['Bid'])

By replicate
------------

.. plot::

    from tbidbaxlipo.plots.nbd_bax_analysis import plot_all_by_replicate
    from tbidbaxlipo.data.parse_bax_bax_fret_nbd_release import df, nbd_residues
    plot_all_by_replicate(df, nbd_residues, ['Release', 'NBD', 'FRET'], activators=['Bid'])

