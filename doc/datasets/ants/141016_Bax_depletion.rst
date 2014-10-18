.. _141016_Bax_depletion:

Bax depletion by liposomes (10/16/14)
=====================================

Preincubation timecourse controls
---------------------------------

.. plot::
    :context:

    from tbidbaxlipo.plots.layout_141016 import \
        plot_release_comparisons, plot_preinc
    plot_preinc()

Timecourse comparison
---------------------

.. plot::
    :context:

    plt.close('all')
    plot_release_comparisons(plot_norm=True, plot_abs=False, bar_plot=True)

