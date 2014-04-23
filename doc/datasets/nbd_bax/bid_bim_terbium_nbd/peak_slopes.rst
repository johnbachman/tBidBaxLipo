Peak slopes
===========

A moving average of the kinetic curves is taken (with a window size of 5) and
then the derivatives are calculated numerically from that.

Bid
---

.. plot::

    from tbidbaxlipo.plots.bid_bim_nbd_release import *
    plot_peak_slopes('Bid')

Bim
---

.. plot::

    from tbidbaxlipo.plots.bid_bim_nbd_release import *
    plot_peak_slopes('Bim')

