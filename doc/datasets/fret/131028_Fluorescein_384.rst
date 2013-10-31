Fluorescein in 384-well plate (10/28/13)
========================================

In this experiment I was interested to see how accurately I could measure
differences in fluorescence intensity at the low nM fluorescein
concentration/signal intensity. The primary motivation was to prepare to do a
competition binding assay with Alexa 488-Bax by FRET. So I took a 384-well
plate (Corning 3575), and in every other row pipetted a 25 nM fluorescein
solution (50 uL). Then I diluted the 25 nM solution by mixing 500 uL of buffer
with 9.5 mL of the solution, and pipetted that 23.75 nM solution into the other
rows of the plate. I then read the plate using the Fluorescein 0.3 sec
protocol on the Wallac, and chose "Costar 384" as the plate layout.

First off, the distribution of means for the 25 nM dye solution wells:

.. plot::
    :context:

    from tbidbaxlipo.plots.layout_131028 import *
    plt.close('all')
    plot(A_timecourses)

Now, the distribution of means for the 23.75 nM dye solution wells:

.. plot::
    :context:

    plt.close('all')
    plot(B_timecourses)

The difference in the peaks of the distributions is apparent, but as we can
see, there are a significant number of outliers for both conditions. We then
look at the edge effects by column and row:

.. plot::
    :context:

    plt.close('all')
    plot_edge_effects(timecourse_wells)

Wow! That is pretty bad. Despite this, the calculated ratio (not even doing
background subtraction) is not too far off, probably because the errors occur
for both the 25 nM and 23.75 nM conditions:

.. ipython:: python

    from tbidbaxlipo.plots.layout_131028 import *
    (A_means, A_sds, A_cvs) = get_means_sds_cvs(A_timecourses)
    A_mean = np.mean(A_means)
    (B_means, B_sds, B_cvs) = get_means_sds_cvs(B_timecourses)
    B_mean = np.mean(B_means)
    print A_mean
    print B_mean
    print B_mean / A_mean

I performed the plate calibration procedure several times, but unfortunately
the edge effects, specifically the drop in fluorescence at the right-most edge,
still persisted. I tried aligning the plate to the upper left corner of the
tray, then the lower right corner; and I looked up the specs for the plate
layout on the Corning website and entered the measurements in. None of these
approaches worked.
