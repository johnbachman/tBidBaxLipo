Bax Activated by Heat, with Auto-Activation, Titration
======================================================

Here we run a stochastic/deterministic comparison of Bax at 43C, with
auto-activation, over the same range of concentrations used for the dataset on
6/14/13. Note that the model uses the one-step approximation for Bax
auto-activation, in which Bax is not saturated:

.. ipython:: python

    from tbidbaxlipo.plots.stoch_det_comparison.bax_heat_auto_titration \
                import jobs, data
    # Macros called by the model
    m = jobs[0].one_cpt_builder().model
    # Rules in the model
    [r.name for r in m.rules]
    # Concentrations in the titration
    bax_concs = [job.params_dict['Bax_0'] for job in jobs]
    bax_concs

We look at a stochastic/deterministic comparison for the minimum and maximum
concentrations. First, for the minimum concentration:

.. plot::
    :context:

    from matplotlib import pyplot as plt
    from tbidbaxlipo.plots.stoch_det_comparison.bax_heat_auto_titration \
                import jobs, data
    from tbidbaxlipo.plots.stoch_det_comparison.plots import \
                plot_bax_timecourse_comparison, plot_dye_release_titration, \
                plot_hist_vs_poisson
    plot_bax_timecourse_comparison(jobs, data, 0)

For the maximum concentration:

.. plot::
    :context:

    plt.close('all')
    plot_bax_timecourse_comparison(jobs, data, len(jobs)-1)

Now, we look at the stochastic/deterministic comparison for the dye release
curves across the titration:

.. plot::
    :context:

    plt.close('all')
    plot_dye_release_titration(jobs, data)

We also look at the distribution of pores across liposomes at several
timepoints:

.. plot::
    :context:

    plt.close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pores', 10)

.. plot::
    :context:

    plt.close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pores', 30)

.. plot::
    :context:

    plt.close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pores', 50)

.. plot::
    :context:

    plt.close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pores', -1)

Two-exponential fits
--------------------

.. plot::
    :context:

    plt.close('all')
    from tbidbaxlipo.plots import bax_heat
    (fmax_arr, k1_arr, k2_arr) = bax_heat.plot_fits_from_CptDataset(data)
    bax_concs = [job.params_dict['Bax_0'] for job in jobs]
    plt.figure()
    plt.plot(bax_concs, fmax_arr)
    plt.figure()
    plt.plot(bax_concs, k1_arr)
    plt.figure()
    plt.plot(bax_concs, k2_arr)


