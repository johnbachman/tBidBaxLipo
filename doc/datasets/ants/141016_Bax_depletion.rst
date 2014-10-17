.. _141016_Bax_depletion:

Bax depletion by liposomes (10/16/14)
=====================================

Preincubation timecourse controls
---------------------------------

.. plot::
    :context:

    from tbidbaxlipo.plots.layout_141016 import \
        lipo_conc_names, plot_release_comparisons, plot_preinc
    plot_preinc()

Timecourse comparison
---------------------

.. plot::
    :context:

    plt.close('all')
    for lipo_conc_name in lipo_conc_names:
        plot_release_comparisons(lipo_conc_name, plot_abs=False)

