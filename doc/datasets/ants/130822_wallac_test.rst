Measuring plate reader error (8/22/13)
======================================

In this experiment I set the Wallac to dispense 20 uL of an ANTS dye solution,
diluted into an appropriate fluorescence range, while taking 0.1 sec reads.
Previously I attempted to do a titration to see if the signal fell into a nice
linear standard curve, but I discovered that above 30-40 uL, pipetting into the
half-well area plates resulted in spilling and wide variance in reads. Hence I
tried 20 uL instead, which seemed to prevent any spilling. I also tried
aligning the plate by its lower-right corner, and this may have helped as well.

I set up the protocol to pipet 20 uL into each well in rows C through H and
take 30 timepoints for each.  The goal was to look at the variability of
measurements just due to signal error (variability within wells)
versus variability due to differences in the amounts of dye pipetted
(variability between wells).

For each timecourse, I truncated the first 10 points, which contained the
signal transient during dispensing. I then took the mean, SD and CV for each
timecourse across the remaining 20 timepoints, and plotted each one in a
histogram. The results show that the **CV due to read error in this range is
generally under 1%**.

.. plot::

    from tbidbaxlipo.plots.layout_130822 import main
    main()

Taking the mean for each well timecourse as the "true" signal for that well,
I took the mean, SD, and CV of these means to look at the pipetting error:

.. ipython:: python

    from tbidbaxlipo.plots.layout_130822 import print_statistics
    print_statistics()

The results show that **the variability between wells due to pipetting error
is considerable, nearly 9%.**

As an additional check, we plot the mean fluorescence values across the rows
to see if there are any edge effects. However, there is no obvious correlation
between the column position and the fluorescence value:

.. plot::

    from tbidbaxlipo.plots.layout_130822 import plot_edge_effects
    plot_edge_effects()
