Timecourses
===========

Raw timecourses from the first experiment, exponential-linear fit
-----------------------------------------------------------------

All of the following plots are for 10uL liposomes.

.. plot::

    from tbidbaxlipo.data.gridv1 import data
    from tbidbaxlipo.grid_analysis import plot_timecourses
    plot_timecourses(data)

Pore data from the first experiment, exponential-linear fit
-----------------------------------------------------------

All of the following plots are for 10uL liposomes.

.. plot::

    from tbidbaxlipo.data.gridv1 import data
    import tbidbaxlipo.grid_analysis as g
    g.plot_timecourses(g.calc_pores(data))

