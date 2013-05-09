Timecourses
===========

Raw timecourses from the first experiment, exponential-linear fit
-----------------------------------------------------------------

All of the following plots are for 10uL liposomes.

.. plot::

    from tbidbaxlipo.data.gridv1 import data
    from tbidbaxlipo.grid_analysis as g
    g.plot_timecourses(data, fittype='biphasic')

Pore data from the first experiment, exponential-linear fit (10uL lipid)
------------------------------------------------------------------------

All of the following plots are for 10uL liposomes.

.. plot::

    from tbidbaxlipo.data.gridv1 import data
    import tbidbaxlipo.grid_analysis as g
    g.plot_timecourses(g.calc_pores(g.calc_bgsub(data)), fittype='biphasic')

Pore data from the second experiment, exponential-linear fit (10uL lipid)
-------------------------------------------------------------------------

All of the following plots are for 10uL liposomes.

.. plot::

    from tbidbaxlipo.data.gridv2 import data
    import tbidbaxlipo.grid_analysis as g
    g.plot_timecourses(g.calc_pores(g.calc_bgsub(data)), fittype='biphasic')

Pore data from the first experiment, exponential-linear fit (3.8uL lipid)
-------------------------------------------------------------------------

All of the following plots are for 3.8uL liposomes.

.. plot::

    from tbidbaxlipo.data.gridv1 import data
    import tbidbaxlipo.grid_analysis as g
    g.plot_timecourses(g.calc_pores(g.calc_bgsub(data)), fixed_conc=3.8, fittype='biphasic')

Pore data from the second experiment, exponential-linear fit (4uL lipid)
-----------------------------------------------------------------------

All of the following plots are for 4uL liposomes.

.. plot::

    from tbidbaxlipo.data.gridv2 import data
    import tbidbaxlipo.grid_analysis as g
    g.plot_timecourses(g.calc_pores(g.calc_bgsub(data)), fixed_conc=4, fittype='biphasic')

Initial rate titration, first experiment
----------------------------------------

.. plot::

    from tbidbaxlipo.data.gridv1 import data
    import tbidbaxlipo.grid_analysis as g
    g.plot_titration(g.calc_initial_slope(g.calc_pores(g.calc_bgsub(data))),'vi')

Initial rate titration, second experiment
-----------------------------------------

.. plot::

    from tbidbaxlipo.data.gridv2 import data
    import tbidbaxlipo.grid_analysis as g
    g.plot_titration(g.calc_initial_slope(g.calc_pores(g.calc_bgsub(data))),'vi')

