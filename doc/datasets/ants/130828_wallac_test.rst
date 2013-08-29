Measuring plate reader error, after repairs (8/28/13)
=====================================================

After noticing the large variability in pipetting in the experiment done
on 8/22/13, I contacted Barb Grant and Kathy Buhl to see if the instrument
could be repaired or calibrated. Kathy Buhl contacted Mike Lemish from
Perkin Elmer, who changed out the syringe and tubing.

Afterwards I ran an experiment very similar to the one on 8/22/13, but this
time I first pipetted 80uL of buffer into every well in the plate, then
programmed the machine to dispense 20 uL into every well prior to reading, then
read the plate 10 times.

Again, I examined the within-well vs. between-well variability.  I took the
mean, SD and CV for each timecourse across the 10 timepoints, and plotted each
one in a histogram. For within-well variability, the results show that the **CV
due to read error in this range is, as before, under 1%**.

.. plot::

    from tbidbaxlipo.plots.layout_130828 import timecourse_wells, plot
    plot(timecourse_wells)

Strangely, the results for between-well variability show that though the
distribution is generally much narrower than before the repairs, there are a
number of "outliers" that have significantly lower readings (notably, there are
none with significantly higher readings).

This has a **strong effect on the CV, making it over 7% despite the tight
distribution among most of the wells:**

.. ipython:: python

    from tbidbaxlipo.plots.layout_130828 import print_statistics, timecourse_wells
    print_statistics(timecourse_wells)

Looking at edge effects, this seems to be due largely to differences in wells
at the left edge of the plate (columns 1-3), though there is one additional
outlier in well D6:

.. plot::

    from tbidbaxlipo.plots.layout_130828 import plot_edge_effects
    plot_edge_effects()

First three columns removed
---------------------------

We replot with the first three columns and well D6 removed. This improves
things considerably:

.. plot::

    from tbidbaxlipo.plots.layout_130828 import row4_wells, plot
    plot(row4_wells)

Removing these outliers **brings the CV down to 1.4%**:

.. ipython:: python

    from tbidbaxlipo.plots.layout_130828 import print_statistics, row4_wells
    print_statistics(row4_wells)

