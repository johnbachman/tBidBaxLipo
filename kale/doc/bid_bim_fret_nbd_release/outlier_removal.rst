Removal of FRET outliers
========================

.. plot::

    from tbidbaxlipo.plots.bid_bim_fret_nbd_release.plot_outliers import \
                plot_fret_outliers, nbd_residues, df, df_pre
    fret_residues = [res for res in nbd_residues if res != 'WT']
    plot_fret_outliers(df, df_pre, ['Bid', 'Bim'], fret_residues)

