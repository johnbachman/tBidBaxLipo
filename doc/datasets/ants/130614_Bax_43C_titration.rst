.. _130614_Bax_43C_titration:

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

.. plot::
    :context:

    plt.close('all')
    plot_fmax_curve(fmax_arr, conc_list)

.. plot::
    :context:

    plt.close('all')
    plot_k1_curve(k1_arr, conc_list)

.. plot::
    :context:

    plt.close('all')
    plot_k2_curve(k2_arr, conc_list)

One-parameter exponential fits
------------------------------

.. plot::

    from matplotlib import pyplot as plt
    from tbidbaxlipo.plots.layout_130614 import df
    from tbidbaxlipo.plots.titration_fits import OneExpNoFmax
    fit = OneExpNoFmax()
    fit.plot_fits_from_dataframe(df)

Linear fits to pore-transformed data
------------------------------------

As expected, produces essentially the same results as the fits to the
one-parameter exponential.

.. plot::

    from matplotlib import pyplot as plt
    from tbidbaxlipo.plots.layout_130614 import pores
    from tbidbaxlipo.util.plate_assay import to_dataframe
    df = to_dataframe(pores)
    from tbidbaxlipo.plots.titration_fits import Linear
    fit = Linear()
    fit.plot_fits_from_dataframe(df)

Two-parameter exponential fits
------------------------------

.. plot::

    from matplotlib import pyplot as plt
    from tbidbaxlipo.plots.layout_130614 import df
    from tbidbaxlipo.plots.titration_fits import OneExpFmax
    fit = OneExpFmax()
    fit.plot_fits_from_dataframe(df)
