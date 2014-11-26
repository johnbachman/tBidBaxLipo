Effect of Bid-membrane saturation on NBD-Bax insertion (11/25/14)
=================================================================

Raw data
--------

Plotting the raw data helps us to identify two outliers, wells E1 and B4.

.. plot::
    :context:

    plt.close('all')
    from tbidbaxlipo.plots.layout_141125 import \
                plot_raw_data, plot_clean_data, plot_bg

    plot_raw_data()

Clean data
----------

Plotting the clean data we see that the outliers have been removed.

.. plot::
    :context:

    plt.close('all')
    plot_clean_data()

Background averages
-------------------

The timecourses of the background wells show that the Bim BH3 actually increases
the background fluorescence in the well somewhat.

.. plot::
    :context:

    plt.close('all')
    plot_bg()
