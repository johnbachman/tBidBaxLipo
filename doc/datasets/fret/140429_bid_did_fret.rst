Bid/DiD liposome FRET, by Justin (4/29/14)
==========================================

1% DiD added to liposomes.

Background controls (BG)
------------------------

Wells containing unlabeled liposomes and a dilution series of unlabeled Bid (no
donor or acceptor). The results show that the unlabeled Bid, even high
concentrations, doesn't substantially contribute to the background. The black
line denotes the average over all of the wells, that will be subtracted from
the donor-only controls.

.. ipython::

    In [1]: from tbidbaxlipo.plots.x140429_Bid_membrane_FRET.calculate_fret import *

    @savefig 140429_bid_did_fret_1.png
    In [7]: bg_avg = get_average_bg(timecourse_wells, bg, plot=True)

Acceptor-only controls (FA)
---------------------------

Wells containing DiD-labeled liposomes, a dilution series of unlabeled Bid, but
not donor. These control for the background of the labeled liposomes in the 568
(donor) channel. The black line denotes the average over all of the wells, that
will be subtracted from the donor + acceptor wells.

.. ipython::

    @savefig 140429_bid_did_fret_2.png
    In [9]: fa_avg = get_average_fa(timecourse_wells, fa, plot=True)

Donor-only (FD)
---------------

Wells containing unlabeled liposomes, a dilution series of unlabeled Bid,
and Bid-568. Used to determine the baseline fluorescence in the absence of FRET.

.. ipython::

    In [2]: fd_well_names = itertools.chain(*fd.values()) #*

    In [3]: fd_wells = extract(fd_well_names, timecourse_wells)

    In [4]: plt.figure()

    @savefig 140429_bid_did_fret_3.png
    In [4]: plot_all(fd_wells)

Here are the average timecourses across replicates:

.. ipython::

    In [1]: (fd_avgs, fd_stds) = averages(timecourse_wells, fd)

    In [2]: plt.figure()

    @savefig 140429_bid_did_fret_4.png
    In [3]: plot_all(fd_avgs)

Donor + acceptor (FDA)
----------------------

Wells containing DiD-labeled liposomes, a dilution series of unlabeled Bid, and Bid-568. Used to measure the FRET, relative to the donor-only condition.

.. ipython::

    In [2]: fda_well_names = itertools.chain(*fda.values()) #*

    In [3]: fda_wells = extract(fda_well_names, timecourse_wells)

    In [4]: plt.figure()

    @savefig 140429_bid_did_fret_5.png
    In [4]: plot_all(fda_wells)

Here are the average timecourses across replicates:

.. ipython::

    In [1]: (fda_avgs, fda_stds) = averages(timecourse_wells, fda)

    In [2]: plt.figure()

    @savefig 140429_bid_did_fret_6.png
    In [3]: plot_all(fda_avgs)

Variability in baseline fluorescence between Bid-568 wells
----------------------------------------------------------

As the timecourses shown above illustrate, the wells with donor (Bid-568), both
the donor-only and donor + acceptor conditions, show substantial differences in
their baseline fluorescence. This cannot be due to increased background in the
presence of high concentrations of unlabeled Bid, since the background plots
(above) also show that the unlabeled Bid doesn't contribute to fluorescence.

If we plot the endpoints of the FD condition as a function of the unlabeled Bid
concentration, we see that the fluorescence goes up and down in an unusual
fashion:

.. ipython::

    @savefig 140429_bid_did_fret_7.png
    In [3]: plot_fd_final_values(timecourse_wells, fd)

I suspect that the discrepancies arose from the a process of doing the dilution
series of the unlabeled Bid into the wells that altered the concentration of
the Bid-568 in a way that was no consistent across wells. It is possible that
this (potential) variability in donor concentration affected the FRET
measurements.


