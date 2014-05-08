Fluorescein read parameters (10/21/13)
======================================

In this experiment I took a 1 uM carboxyfluorescein stock solution and did
three two-fold serial dilutions in a plate (Corning 3686). The goal was to see
if the serial dilutions using the Chaney adapter and my new jig for holding the
syringe would improve reproducibility in the serial dilution.  In addition, I
wanted to see how plate reading parameters would affect the error at different
fluorescein concentrations/signal intensities.

After doing the dilutions (50 uL in each well) I added 50 uL of buffer into
each well to mimic the conditions of a Bax/Bid FRET experiment, which reduced
the concentrations by half.

I took 20 reads of the whole plate using the **Fluorescein 1 sec standard
protocol.** To see data from other read parameters (0.1 sec or 0.3 sec), edit
tbidbaxlipo.plots.layout_131021.py to load the appropriate data.  The 0.3 sec
custom protocol was otherwise unchanged from the standard protocol.

First, we look at the data over time, both separately for each dilution series
and also averaged. The readings are fairly stable over time and the error bars
between the replicates are small.

.. plot::
    :context:

    plt.close('all')
    from tbidbaxlipo.plots.layout_131021 import *
    plot_data()

Next, we take the data and for each well, take the mean and SD for the 20 reads
and plot as a bar plot. Viewed this way, the relative error can be seen to be quite small.

.. plot::
    :context:

    plt.close('all')
    plot_bar_plot()

To examine this more closely, we plot the CV as a function of concentration and
can see that it rapidly decays as the concentration/signal increases, down
to a level of about 0.4 for high concentrations:

.. plot::
    :context:

    plt.close('all')
    plot_cvs()

Now we plot the dilution series one top of each other to look at the
reproducibility of the different dilution series. It can be seen here that the
dilution process is highly reproducible:

.. plot::
    :context:

    plt.close('all')
    plot_dilution_series()

Finally, we fit the fluorescence titration to a variety of curves, and find that
the data is best fit not by a line, but rather by a power law:

.. plot::
    :context:

    plt.close('all')
    plot_fits()

