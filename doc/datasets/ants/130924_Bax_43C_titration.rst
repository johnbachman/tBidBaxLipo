Bax Titration at 43C (9/24/13)
==============================

In this experiment I used more liposomes (2 mL + 500 uL buffer, for an
80% solution) and less protein (dilution series down from 126 nM Bax.

Raw data
--------

.. plot::

    from tbidbaxlipo.plots.layout_130924 import plot_data
    plot_data()

Two-exponential fits
--------------------

Fits of each replicate to the three-parameter equation

.. math::

    F_{max} \left(1 - e^{-k_1 (1 - e^{-k_2 t})t}\right)

.. plot::
    :context:

    from tbidbaxlipo.plots.layout_130924 import *
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

    from tbidbaxlipo.plots.layout_130924 import df
    from tbidbaxlipo.plots.titration_fits import OneExpNoFmax
    fit = OneExpNoFmax()
    fit.plot_fits_from_dataframe(df)

Linear fits to pore-transformed data
------------------------------------

.. plot::

    from tbidbaxlipo.plots.layout_130924 import pores
    from tbidbaxlipo.util.plate_assay import to_dataframe
    df = to_dataframe(pores)
    from tbidbaxlipo.plots.titration_fits import Linear
    fit = Linear()
    fit.plot_fits_from_dataframe(df)

Two-parameter exponential fits
------------------------------

.. plot::

    from tbidbaxlipo.plots.layout_130924 import df
    from tbidbaxlipo.plots.titration_fits import OneExpFmax
    fit = OneExpFmax()
    fit.plot_fits_from_dataframe(df)

Two-exponential fits
--------------------

.. plot::

    from matplotlib import pyplot as plt
    from tbidbaxlipo.plots.layout_130924 import df
    from tbidbaxlipo.plots.titration_fits import TwoExp
    fit = TwoExp()
    fit.plot_fits_from_dataframe(df)

Two-exponential fits, active subset
-----------------------------------

.. plot::

    from matplotlib import pyplot as plt
    from tbidbaxlipo.plots.layout_130924 import subset_df
    from tbidbaxlipo.plots.titration_fits import TwoExp
    fit = TwoExp()
    fit.plot_fits_from_dataframe(subset_df)

