Plots and analysis of NBD data
==============================

Raw data
--------

.. plot::

    from tbidbaxlipo.nbd_analysis import *
    plot_raw()

Normalized to [0, 1]
--------------------

.. plot::

    from tbidbaxlipo.nbd_analysis import *
    plot_normalized()

Averaged across replicates
--------------------------

.. plot::

    from tbidbaxlipo.nbd_analysis import *
    plot_avg(normalize=False)

Normalized and averaged across replicates
-----------------------------------------

.. plot::

    from tbidbaxlipo.nbd_analysis import *
    plot_avg(normalize=True)

Fits with mathematical functions
--------------------------------

.. plot::

    from tbidbaxlipo.nbd_analysis import *
    plot_fit()

