Initial rates
=============

The first n points (in this case 4) are fitted with a line to determine
the initial slope. The error bars are the standard error of the fit
parameters returned by the fitting function. When the ratio is taken,
the standard error of the ratio is estimated by sampling.

.. plot::

    from tbidbaxlipo.plots.bid_bim_nbd_release import *
    plot_initial_rates('Bid')
    plot_initial_rates('Bim')


