Fluorescein dilution series on FlexStation (10/30/13)
=====================================================

To get to know the capabilities of the FlexStation 3 better I set up a
fluorescein dilution series with 5 replicates. The idea was to optimize read
settings and signal-to-noise for concentrations in the range of 1-10 nM
fluorescein.

I prepared a 200 nM stock by adding 2 uL of 100 uM stock to 998 uL of buffer.
I did the dilution in a half-well area black plate (Corning 3686) starting with
100 uL in the first well and doing two-fold dilutions. I did five parallel
dilutions in rows A-E.

For the read I tried three settings:

* medium sensitivity, 6 reads (the default)
* high sensitivity, 22 reads (this didn't increase read time too much, but
  resulted in saturation of the 100 nM concentration)
* medium sensitivity, 100 reads

The main conclusions are:

* The fluorescein signal is detectable and linear down to concentrations of
  100 picomolar.
* Increased PMT sensitivity decreases CV of reads at low sensitivity if
  saturation is avoided.
* A moderate increase in number of reads can reduce CV without severely slowing
  down the read time.

Medium sensitivity, 6 reads
---------------------------

First, the raw data, showing pretty good stability with a handful of outliers:

.. plot::
    :context:

    from tbidbaxlipo.plots.layout_131030 import *
    close('all')
    [timecourse_wells, timecourse_averages, timecourse_stds,
     conc_list, means, stds, log_concs, log_means, log_stds] = \
                    get_wells('131030_Fluorescein_Flex_Med_6read.txt', layout)
    plot_data(timecourse_wells, timecourse_averages, timecourse_stds)

The bar plot shows that the within-read well variability is small relative
to the between-well variability; there is also a striking increase in signal
in going from row A to row E. Is this the result of pipetting error (I didn't
change the tips) or is there some kind of edge effect in the plate?

.. plot::
    :context:

    close('all')
    plot_bar_plot(means, stds)

The CVs are fairly high within-well CV at the low concentrations, above 20%
even!

.. plot::
    :context:

    close('all')
    plot_cvs(means, stds, log_concs)

The dilution series appears nicely linear:

.. plot::
    :context:

    close('all')
    plot_dilution_series(conc_list, means, stds, log_concs, log_means, log_stds)

The linear fit shows this as well:

.. plot::
    :context:

    close('all')
    plot_fits(means, log_means, conc_list, log_concs)

The correlation coefficient is very good:

.. ipython:: python

    from tbidbaxlipo.plots.layout_131030 import *
    [timecourse_wells, timecourse_averages, timecourse_stds,
     conc_list, means, stds, log_concs, log_means, log_stds] = \
                    get_wells('131030_Fluorescein_Flex_Med_6read.txt', layout)
    plot_fits(means, log_means, conc_list, log_concs)

High sensitivity, 22 reads
--------------------------

The CVs were good for these read parameters but the 100 nM concentration was
saturated.

.. plot::
    :context:

    close('all')
    [timecourse_wells, timecourse_averages, timecourse_stds,
     conc_list, means, stds, log_concs, log_means, log_stds] = \
                    get_wells('131030_Fluorescein_Flex_High_22read.txt', layout)
    plot_data(timecourse_wells, timecourse_averages, timecourse_stds)

.. plot::
    :context:

    close('all')
    plot_bar_plot(means, stds)

The CVs appeared to be the lowest for these read parameters:

.. plot::
    :context:

    close('all')
    plot_cvs(means, stds, log_concs)

.. plot::
    :context:

    close('all')
    plot_dilution_series(conc_list, means, stds, log_concs, log_means, log_stds)

.. plot::
    :context:

    close('all')
    plot_fits(means, log_means, conc_list, log_concs)

.. ipython:: python

    [timecourse_wells, timecourse_averages, timecourse_stds,
     conc_list, means, stds, log_concs, log_means, log_stds] = \
                    get_wells('131030_Fluorescein_Flex_High_22read.txt', layout)
    plot_fits(means, log_means, conc_list, log_concs)

Medium sensitivity, 100 reads
-----------------------------

The additional reads significantly increased reading time, with no proportional
improvement in the CVs at low concentrations.

.. plot::
    :context:

    close('all')
    [timecourse_wells, timecourse_averages, timecourse_stds,
     conc_list, means, stds, log_concs, log_means, log_stds] = \
                    get_wells('131030_Fluorescein_Flex_Med_100read.txt', layout)
    plot_data(timecourse_wells, timecourse_averages, timecourse_stds)

.. plot::
    :context:

    close('all')
    plot_bar_plot(means, stds)

.. plot::
    :context:

    close('all')
    plot_cvs(means, stds, log_concs)

.. plot::
    :context:

    close('all')
    plot_dilution_series(conc_list, means, stds, log_concs, log_means, log_stds)

.. plot::
    :context:

    close('all')
    plot_fits(means, log_means, conc_list, log_concs)

.. ipython:: python

    [timecourse_wells, timecourse_averages, timecourse_stds,
     conc_list, means, stds, log_concs, log_means, log_stds] = \
                    get_wells('131030_Fluorescein_Flex_Med_100read.txt', layout)
    plot_fits(means, log_means, conc_list, log_concs)

