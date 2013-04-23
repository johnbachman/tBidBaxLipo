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

Single exponential model
~~~~~~~~~~~~~~~~~~~~~~~~

Fits to the equation

.. math::

    F_0 + F_{max}\left(1 - e^{-kt}\right)

.. plot::

    from tbidbaxlipo.nbd_analysis import *
    plot_fit(fittype='single_exp')


Exponential-linear model
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    F_0 + F_{max}\left(1 - e^{-kt}\right) + mt

.. plot::

    from tbidbaxlipo.nbd_analysis import *
    plot_fit(fittype='exp_lin')

Two-exponential model
~~~~~~~~~~~~~~~~~~~~~

.. math::

    F_0 + F_{max_1}\left(1 - e^{-k_1 t}\right) +
    F_{max_2}\left(1 - e^{-k_2 t}\right)

.. plot::

    from tbidbaxlipo.nbd_analysis import *
    plot_fit(fittype='double_exp')


