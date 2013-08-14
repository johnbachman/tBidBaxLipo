Bax Activated by Heat
=====================

.. plot::
    :context:

    from tbidbaxlipo.plots.stoch_det_comparison.bax_heat import jobs, data
    from tbidbaxlipo.plots.stoch_det_comparison.plots import *
    plot_bax_timecourse_comparison(jobs, data, 0)

.. plot::
    :context:

    close('all')
    plot_dye_release_titration(jobs, data)

.. plot::
    :context:

    close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pores', -1)

