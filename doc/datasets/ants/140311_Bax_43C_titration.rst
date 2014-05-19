.. _140311_Bax_43C_titration:

Bax Titration at 43C (3/11/14)
==============================

Replicate of 6/14/13 **TODO** add link here.

Raw data
--------

.. plot::
    :context:

    from tbidbaxlipo.plots.layout_140311 import *
    plt.close('all')
    plot_data()

Initial rates
-------------

To see how the initial rate of permeabilization scales with the Bax
concentration, we fit the first 50 points to a line (through the origin), and
calculate the slope, then plot the fitted slopes as a function of Bax
concentration. First, the results from the fits to the background subtracted,
normalized data:

.. plot::
    :context:

    plt.close('all')
    (k1_arr, conc_list) = plot_linear_fits(plot=True, max_time=50)

Now we plot the initial slopes as a function of the concentration, and attempt
to fit this rate plot with an exponential curve or a linear fit to the first
several points:

.. plot::
    :context:

    plt.close('all')
    bax_rate_slope = plot_initial_slope_curves(k1_arr, conc_list)

Two-exponential fits
--------------------

Now we fit with the two-exponential equation from Kushnareva et al.:

.. math::

    F_{max} \left(1 - e^{-k_1 (1 - e^{-k_2 t})t}\right)

.. plot::
    :context:

    plt.close('all')
    (fmax_arr, k1_arr, k2_arr, conc_list) = \
           titration_fits.plot_two_exp_fits(bgsub_norm_wells, layout, plot=True)
    titration_fits.plot_fmax_curve(fmax_arr, conc_list)
    titration_fits.plot_k1_curve(k1_arr[:,1:], conc_list[1:])
    titration_fits.plot_k2_curve(k2_arr[:,1:], conc_list[1:])


