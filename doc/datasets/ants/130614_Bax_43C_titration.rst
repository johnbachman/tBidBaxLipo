Bax Titration at 43C (6/14/13)
==============================

The idea behind this experiment was to look at the dynamics of Bax
permeabilization in a unimolecular scenario (with no tBid) over a range of
concentrations to create a dataset that could be used to constrain possible
mechanisms.

Raw data
--------

.. plot::

    from tbidbaxlipo.plots.layout_130614 import plot_data
    plot_data()

Two-exponential fits
--------------------

Fits of each replicate to the three-parameter equation

.. math::

    F_{max} \left(1 - e^{-k_1 (1 - e^{-k_2 t})t}\right)

.. plot::
    :context:

    from tbidbaxlipo.plots.layout_130614 import *
    (fmax_arr, k1_arr, k2_arr, conc_list) = plot_two_exp_fits()

:math:`F_{max}` curves
----------------------

.. plot::
    :context:

    plt.close('all')
    plot_fmax_curve(fmax_arr, conc_list)

:math:`k_1` curves
------------------

.. plot::
    :context:

    plt.close('all')
    plot_k1_curve(k1_arr, conc_list)

:math:`k_2` curves
------------------

.. plot::
    :context:

    plt.close('all')
    plot_k2_curve(k2_arr, conc_list)

