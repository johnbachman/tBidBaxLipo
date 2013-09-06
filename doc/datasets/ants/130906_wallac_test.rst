Measuring plate reader error, full-area plate (9/6/13)
======================================================

The pipetting experiment done on :ref:`8/22/13 <130822_wallac_test>` revealed that
pipetting more than 20uL or so into a half-well area plate (Corning 3686) led
to spillage and significant errors.

On the other hand, the pipetting experiment done on :ref:`8/28/13
<130828_wallac_test>` showed significant edge effects, presumably due to error
in pipetting at the left-edge of the plate (though I did not rule out errors
due to the position of the fluorescence read).  I expected that this was due in
part to the half-well area plate leading to missed pipetting or other errors.
With this in mind, I purchased full-area 96-well area plates (with the
non-binding surface coating, Corning 3991) and tested these plates for
pipetting reproducibility.

Anticipating that the full-well area plates might not have the spillage problem
with volumes of 50 uL or so, I pipetted 50 uL of buffer into each well and then
had the Wallac pipette 50 uL of an ANTS dye solution into each well, then
read the plate 20 times.

Again, I examined the within-well vs. between-well variability.  I took the
mean, SD and CV for each timecourse across the 10 timepoints, and plotted each
one in a histogram.

The results show the following:

* The very first well, A01, is significantly lower than the others. I have seen
  this in other pipetting experiments and expect that this has to due
  essentially with priming the syringe with fluid.  The solution is probably to
  just avoid making measurements from the first well.
* There is a time dependence in the fluorescence reads, as they appear to go up
  over the first half before plateauing (or even going down).
* Even including this variability due to the increase/decrease, the within-well
  variability is under 1% CV on average.

.. plot::

    from tbidbaxlipo.plots.layout_130906 import timecourse_wells, plot
    plot(timecourse_wells)

Encouragingly, the results for between-well variability show that even though
there are a tiny handful of outliers (well A01, for example) **the pipetting for
50 uL into this plate is highly reproducible, with a CV even less than that of
the within-well variability, of 0.79%!**

.. ipython:: python

    from tbidbaxlipo.plots.layout_130906 import print_statistics, timecourse_wells
    print_statistics(timecourse_wells)

There do appear to be edge effects, however--even leaving out the outlier A01,
there is a trend of low-high-low across the plate. This should probably be addressed through the plate calibration procedure.

.. plot::

    from tbidbaxlipo.plots.layout_130906 import plot_edge_effects
    plot_edge_effects()


