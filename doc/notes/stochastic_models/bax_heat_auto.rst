Bax Activated by Heat, with Auto-Activation
===========================================

.. ipython:: python

    from tbidbaxlipo.plots.stoch_det_comparison.bax_heat_auto import jobs, data
    m = jobs[0].one_cpt_builder().model
    [r.name for r in m.rules]

.. plot::
    :context:

    from tbidbaxlipo.plots.stoch_det_comparison.bax_heat_auto import jobs, data
    from tbidbaxlipo.plots.stoch_det_comparison.plots import *
    plot_bax_timecourse_comparison(jobs, data, 0)

.. plot::
    :context:

    close('all')
    plot_dye_release_titration(jobs, data)

.. plot::
    :context:

    close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pores', 10)

.. plot::
    :context:

    close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pores', 30)

.. plot::
    :context:

    close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pores', 50)

.. plot::
    :context:

    close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pores', -1)

