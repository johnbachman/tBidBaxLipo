NBD-Bax Kinetics, PTI
=====================

An NBD-Bax dataset provided by David Andrews in the fall of 2011.  Includes
kinetics of several NBD-Bax mutants measured by the PTI fluorimeter, with high
precision readings every two seconds.

Includes data for the following single-cysteine mutants (three replicates
each):

* c3
* c62
* c120
* c122
* c126

In this dataset, Bax was first equilibrated with liposomes to determine the
background fluorescence, then cBid was subsequently added for measurement of
the insertion kinetics. However, the Andrews lab later determined that there is
a lag in cBid activity due to a slow conformational change in Bid
[ShamasDin2013]_, so the curves in this dataset may manifest delayed activation
as a result.

Raw data
--------

The "raw data" here are curves in which the fluorescence increase has been
normalized to the baseline (pre-cBid), so the y-axis values represent the
fold-change increase.

.. plot::

    from tbidbaxlipo.plots.nbd_analysis import *
    plot_raw()

Normalized to [0, 1]
--------------------

The problem with simply normalizing each curve to [0, 1] based on the minimum
and maximum values of each is that outlier values (spikes in the data) create
an artificially high maximum value and distort the normalization. To get around
this, here I chose one curve to normalize to [0, 1] based on its min and max
values, then chose fit a scaling parameter and an offset to "stretch" the other
curves to it via a linear transformation. The fact that the curves overlap very
closely when normalized in this way shows that the kinetics are remarkably
reproducible, even if the magnitudes (and fold-change magnitudes) are not.

.. plot::

    from tbidbaxlipo.plots.nbd_analysis import *
    plot_normalized()

Averaged across replicates
--------------------------

.. plot::

    from tbidbaxlipo.plots.nbd_analysis import *
    plot_avg(normalize=False)

Normalized and averaged across replicates
-----------------------------------------

.. plot::

    from tbidbaxlipo.plots.nbd_analysis import *
    plot_avg(normalize=True)

Fits with mathematical functions
--------------------------------

Single exponential model
~~~~~~~~~~~~~~~~~~~~~~~~

Fits to the equation

.. math::

    F_0 + F_{max}\left(1 - e^{-kt}\right)

.. plot::

    from tbidbaxlipo.plots.nbd_analysis import *
    plot_fit(fittype='single_exp')


Exponential-linear model
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    F_0 + F_{max}\left(1 - e^{-kt}\right) + mt

.. plot::

    from tbidbaxlipo.plots.nbd_analysis import *
    plot_fit(fittype='exp_lin')

Two-exponential model
~~~~~~~~~~~~~~~~~~~~~

.. math::

    F_0 + F_{max_1}\left(1 - e^{-k_1 t}\right) +
    F_{max_2}\left(1 - e^{-k_2 t}\right)

.. plot::

    from tbidbaxlipo.plots.nbd_analysis import *
    plot_fit(fittype='double_exp')


