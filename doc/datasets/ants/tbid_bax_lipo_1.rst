ANTS Release, tBid/Bax/Lipo Titration, v1
=========================================

Raw data, exponential-linear fit
--------------------------------

All of the following plots are for 10uL liposomes.

.. plot::

    from tbidbaxlipo.data.gridv1 import data
    import tbidbaxlipo.plots.grid_analysis as g
    g.plot_timecourses(data, fittype='biphasic')

Pore-transformed data, exponential-linear fit (10uL lipid)
----------------------------------------------------------

All of the following plots are for 10uL liposomes.

.. plot::

    from tbidbaxlipo.data.gridv1 import data
    import tbidbaxlipo.plots.grid_analysis as g
    g.plot_timecourses(g.calc_pores(g.calc_bgsub(data)), fittype='biphasic')

Pore-transformed data, exponential-linear fit (3.8uL lipid)
-----------------------------------------------------------

All of the following plots are for 3.8uL liposomes.

.. plot::

    from tbidbaxlipo.data.gridv1 import data
    import tbidbaxlipo.plots.grid_analysis as g
    g.plot_timecourses(g.calc_pores(g.calc_bgsub(data)), fixed_conc=3.8, fittype='biphasic')

Initial rate titration
----------------------

.. plot::

    from tbidbaxlipo.data.gridv1 import data
    import tbidbaxlipo.plots.grid_analysis as g
    g.plot_titration(g.calc_initial_slope(g.calc_pores(g.calc_bgsub(data))),'vi')
