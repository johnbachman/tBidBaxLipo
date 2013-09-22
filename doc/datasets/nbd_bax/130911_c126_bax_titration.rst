.. _130911_c126_bax_titration:

NBD-c126-Bax titration (9/11/13)
================================

Data from Justin Kale, dated 9/11/13. His comments:

    Hi John,

    Attached is the data for the Bax 126C-NBD titration. I used 30 ul of 2mg/ml
    liposomes for a final concentration of 0.6 mg/ml liposomes or 9.275 nM. Bax
    was titrated from 800-12.5 nM, however looking at the fluorescence data the
    12.5 seems to have some pipetting error as fluorescence counts are similar
    to 25 nM. I attempted to pipette 0.1 uL so this is likely why.

    Also, the software crashed between the background read and the fluorescence
    read once Bax was added. So there are no background reads for each well (I
    used the wells with 0 Bax as backgrounds for every well). Additionally it
    seems a good chunk of the initial curves for the lower Bax concentrations
    were missed, as evident by the higher F/F0 values for the initial point.

    For the wells with cBid I pre-incubated the liposomes with cBid at 2 mg/ml,
    and then added 30 ul of these liposomes to each well for the final
    concentration of 0.6 mg/ml.

    As suspected higher Bax concentrations slows down the kinetics of
    insertion, and it appears they all reach a similar F/F0 (by just looking at
    the data)

    If you have any further questions, let me know.

    Thanks,

    Justin

:math:`F_0` values
------------------

First we look at the initial values for both the 0 nM and 20 nM cBid cases.
We fit the :math:`F_0` values to a line to see if the initial fluorescence
scales linearly with the amount of NBD-Bax added.

For the no-cBid case, the initial values lie almost perfectly on a line, with
an :math:`R^2` of 0.9995. For the 20 nM cBid case, the values are linear,
**except that they do not pass through the point for 0 Bax.**  The fit shown
in the plot below excludes the 0 Bax point when fitting the 20 nM cBid case.
This suggests that an initial portion of the timecourses for 20 nM Bax was
missed, but that the increase in fluorescence during the missed period was
itself linear with Bax concentration.

.. ipython::

    In [1]: from tbidbaxlipo.data.parse_09_11_2013_Bax_126_NBD_titration \
       ...: import *

    In [2]: plot_f0()

    @suppress
    In [3]: plt.savefig('_static/130911_c126_bax_titration_1.png')

.. image:: ../../_static/130911_c126_bax_titration_1.png
    :width: 6in

Raw fluorescence timecourses
----------------------------

For reference, we plot the timecourses for both the 0 and 20 nM cBid
conditions:

.. ipython::

    In [4]: plot_raw_timecourses(bid_conc=0)

    @suppress
    In [5]: plt.savefig('_static/130911_c126_bax_titration_2.png')

.. image:: ../../_static/130911_c126_bax_titration_2.png
    :width: 6in

.. ipython::

    In [6]: plot_raw_timecourses(bid_conc=20)

    @suppress
    In [7]: plt.savefig('_static/130911_c126_bax_titration_3.png')

.. image:: ../../_static/130911_c126_bax_titration_3.png
    :width: 6in

Background subtracted, normalized by 0 cBid :math:`F_0`
-------------------------------------------------------

Now we look at the data after subtracting out the baseline fluorescence
increase in the no Bax (liposomes and Bid only) condition, normalized by the
initial fluorescence values from the 0 cBid timecourses.

Though the no-Bid timecourses normalized in this way are obviously noisy due
to the low signal, they appear to mostly lie on top of each other, suggesting
that the kinetics of spontaneous insertion scale uniformly with concentration,
with no saturation. It will be interesting to see if this holds true in
an experiment with heated Bax.

.. ipython::

    In [8]: plot_normalized(bid_conc=0, bid_conc_for_normalization=0, \
       ...: subtract_background=True)

    @suppress
    In [9]: plt.savefig('_static/130911_c126_bax_titration_4.png')

.. image:: ../../_static/130911_c126_bax_titration_4.png
    :width: 6in

With cBid added, from the shape of the curves it appears that there is a large
portion of the initial timecourses missing. That said, it appears that the
steady-state F/F0 values are fairly constant across concentrations, with the
exception of the 12.5 nM condition, which is lower than the others, and the 25
nM condition which is higher (though both have a lot of error).  Moreover, it
appears that the kinetics slow down at higher Bax concentrations, suggestive of
a saturation effect.

.. ipython::

    In [8]: plot_normalized(bid_conc=20, bid_conc_for_normalization=0, \
       ...: subtract_background=True)

    @suppress
    In [9]: plt.savefig('_static/130911_c126_bax_titration_5.png')

.. image:: ../../_static/130911_c126_bax_titration_5.png
    :width: 6in

To evaluate this accurately, we fit the curves with parameters describing
the rate and maximum steady state insertion and plot the values of these
parameters as a function of Bax concentration. However, to fit the curves, we
have to account for the lag before the start of the measurement, which we can
do by fitting an additional parameter :math:`t_0` describing the estimated
length of the lag. We use the equation :math:`1 + F_{max} (1 - e^{-k_1 (t +
t_0)})` and get the following fits:

.. ipython::

    In [8]: (k1, fmax, t0) = plot_normalized(bid_conc=20,
       ...: bid_conc_for_normalization=0, subtract_background=True,
       ...: do_fit=True, t0_val=None)

    @suppress
    In [9]: plt.savefig('_static/130911_c126_bax_titration_6.png')

.. image:: ../../_static/130911_c126_bax_titration_6.png
    :width: 6in

The values for the fitted parameters are as follows:

.. ipython:: python

    # Format parameter values in a table
    tt = Texttable()
    tt.header(['[Bax]', 'k1', 'Fmax', 't0'])
    tt.set_cols_dtype(['f', 'e', 'f', 'f'])
    tt.add_rows(reversed(zip(bax_concs, k1, fmax, t0)), header=False)
    print tt.draw()

The fitted values for the delay parameter :math:`t_0` show a pattern, with the
higher Bax concentrations better fit by shorter delays, and the lower
concentrations with longer delays (note also that the poor fit to 0 Bax is
irrelevant). This appears to be an artifact of the shape of the high
concentration curves, as Justin actually pipetted the Bax in the reverse order
(highest concentrations first).

Turning to the parameter plots, the plot for :math:`k_1` vs. Bax shows a very
clear slowing down in insertion rate as the Bax concentration is increased:

.. ipython::

    In [1]: plt.figure()

    In [2]: plt.plot(bax_concs[:-1], k1[:-1], linewidth=2, marker='o', color='r')

    @suppress
    In [3]: plt.xlabel('[Bax] (nM)')

    @suppress
    In [4]: plt.ylabel('$k_1$')

    @suppress
    In [5]: plt.title('$k_1$ vs. [Bax]')

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_7.png')

.. image:: ../../_static/130911_c126_bax_titration_7.png
    :width: 6in

On the other hand, the :math:`F_{max}` values, though a bit noisy, hover
generally in the range of 3.0-3.3 without an obvious pattern. It would be
a good idea to get error bars on these before making any strong conclusions.

.. ipython::

    In [1]: plt.figure()

    In [2]: plt.plot(bax_concs[:-1], fmax[:-1], linewidth=2, marker='o', color='r')

    @suppress
    In [3]: plt.xlabel('[Bax] (nM)')

    @suppress
    In [4]: plt.ylabel('$F_{max}$')

    @suppress
    In [5]: plt.title('$F_{max}$ vs. [Bax]')

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_8.png')

.. image:: ../../_static/130911_c126_bax_titration_8.png
    :width: 6in

Background-subtracted, normalized by 0 nM cBid, fit with fixed :math:`t_0`
--------------------------------------------------------------------------

A potential concern about the results shown above is that the fitted values for
the parameters :math:`k_1` and :math:`F_{max}` may be influenced by the fitted
value for the parameter :math:`t_0`, which is allowed to vary for each curve.
Some of the fitted delays seem unrealistically long (e.g. ~700-800 sec for 12.5
and 25 nM Bax), and it also seems likely that a uniform delay for all wells due
to the software crash that Justin described would be large relative to the
smaller variable delays in pipetting the different wells.

To address this, we re-fit the data (normalized by the 0 nM cBid :math:`F_0`
values) while fixing :math:`t_0` to a constant value. We fix it to the average
over all the fitted values found in the previous step (ignoring the
fitted value for 0 Bax):

.. ipython::

    In [1]: mean_t0 = np.mean(t0[:-1]); print mean_t0

    In [8]: (k1, fmax, t0) = plot_normalized(bid_conc=20,
       ...: bid_conc_for_normalization=0, subtract_background=True,
       ...: do_fit=True, t0_val=mean_t0)

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_9.png')

.. image:: ../../_static/130911_c126_bax_titration_9.png
    :width: 6in

Fortunately, the fits still appear to be quite good despite the fact that we
have removed a parameter for each curve. The initial ramp-up for the higher
concentration curves is less well fit but this is typical for single
exponential fits of these curves (and as such is actually a bit reassuring).
The re-fit values for :math:`k_1` and :math:`F_{max}` are as follows:

.. ipython:: python

    # Format parameter values in a table
    tt = Texttable()
    tt.header(['[Bax]', 'k1', 'Fmax', 't0'])
    tt.set_cols_dtype(['f', 'e', 'f', 'f'])
    tt.add_rows(reversed(zip(bax_concs, k1, fmax, t0)), header=False)
    print tt.draw()

The titration plot for :math:`k_1` still shows a very clear reduction in the
insertion rate with increasing Bax; in fact some of the noise appears to
have been reduced:

.. ipython::

    In [1]: plt.figure()

    In [2]: plt.plot(bax_concs[:-1], k1[:-1], linewidth=2, marker='o', color='r')

    @suppress
    In [3]: plt.xlabel('[Bax] (nM)')

    @suppress
    In [4]: plt.ylabel('$k_1$')

    @suppress
    In [5]: plt.title('$k_1$ vs. [Bax]')

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_10.png')

.. image:: ../../_static/130911_c126_bax_titration_10.png
    :width: 6in

The :math:`F_{max}` values unchanged, which makes sense since the steady-state
value shouldn't be much affected by the slight change in start time:

.. ipython::

    In [1]: plt.figure()

    In [2]: plt.plot(bax_concs[:-1], fmax[:-1], linewidth=2, marker='o', color='r')

    @suppress
    In [3]: plt.xlabel('[Bax] (nM)')

    @suppress
    In [4]: plt.ylabel('$F_{max}$')

    @suppress
    In [5]: plt.title('$F_{max}$ vs. [Bax]')

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_11.png')

.. image:: ../../_static/130911_c126_bax_titration_11.png
    :width: 6in

Normalized by 0 cBid :math:`F_0`, no background subtraction
-----------------------------------------------------------

The following plots demonstrate the importance of background subtraction in
plotting and fitting the curves correctly. Here we normalize the raw values by
the corresponding initial fluorescence :math:`F_0` from the no Bid titration,
but we do not subtract the background (no Bax) fluorescence values for each
timepoint before doing so.

In the plot of the no-Bid timecourses we can clearly see that the 0 Bax
condition, once normalized, has a non-negligible increase in baseline
fluorescence. In addition, it appears that the kinetics of the other curves
scale with concentration, with the lower concentrations increasing more slowly
than the higher concentrations. However, this is merely due to the fact that
since we have not subtracted the background from the :math:`F_0` values, the
fold-change increase over background appears to be less for the curves with low
signal (low concentration).

.. ipython::

    In [8]: plot_normalized(bid_conc=0, bid_conc_for_normalization=0, \
       ...: subtract_background=False)

    @suppress
    In [9]: plt.savefig('_static/130911_c126_bax_titration_12.png')

.. image:: ../../_static/130911_c126_bax_titration_12.png
    :width: 6in

Similarly, the plots for the 20 nM cBid condition seem to show that the steady
state :math:`F/F_0` values go up in a saturating fashion with concentration:

.. ipython::

    In [8]: plot_normalized(bid_conc=20, bid_conc_for_normalization=0, \
       ...: subtract_background=False)

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_13.png')

.. image:: ../../_static/130911_c126_bax_titration_13.png
    :width: 6in

The fits for the :math:`F_{max}` parameter show the apparent increase in
steady-state fluorescence, which looks deceptively like a binding curve:

.. ipython::

    In [8]: (k1, fmax, t0) = plot_normalized(bid_conc=20,
       ...: bid_conc_for_normalization=0, subtract_background=False,
       ...: do_fit=True, t0_val=None)

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_14.png')

.. image:: ../../_static/130911_c126_bax_titration_14.png
    :width: 6in

.. ipython:: python

    # Format parameter values in a table
    tt = Texttable()
    tt.header(['[Bax]', 'k1', 'Fmax', 't0'])
    tt.set_cols_dtype(['f', 'e', 'f', 'f'])
    tt.add_rows(reversed(zip(bax_concs, k1, fmax, t0)), header=False)
    print tt.draw()

.. ipython::

    In [1]: plt.figure()

    In [2]: plt.plot(bax_concs[:-1], k1[:-1], linewidth=2, marker='o', color='r')

    @suppress
    In [3]: plt.xlabel('[Bax] (nM)')

    @suppress
    In [4]: plt.ylabel('$k_1$')

    @suppress
    In [5]: plt.title('$k_1$ vs. [Bax]')

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_15.png')

.. image:: ../../_static/130911_c126_bax_titration_15.png
    :width: 6in

.. ipython::

    In [1]: plt.figure()

    In [2]: plt.plot(bax_concs[:-1], fmax[:-1], linewidth=2, marker='o', color='r')

    @suppress
    In [3]: plt.xlabel('[Bax] (nM)')

    @suppress
    In [4]: plt.ylabel('$F_{max}$')

    @suppress
    In [5]: plt.title('$F_{max}$ vs. [Bax]')

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_16.png')

.. image:: ../../_static/130911_c126_bax_titration_16.png
    :width: 6in

Normalized by 20 nM cBid :math:`F_0`
------------------------------------

Alternatively, we can normalize the 20 nM cBid timecourses by their own initial value, which produces the following curves as a result:

.. ipython::

    In [8]: plot_normalized(bid_conc=20, bid_conc_for_normalization=20, \
       ...: subtract_background=True)

    @suppress
    In [13]: plt.savefig('_static/130911_c126_bax_titration_17.png')

.. image:: ../../_static/130911_c126_bax_titration_17.png
    :width: 6in

Despite the fact that we've subtracted the fluorescence from the no-Bax
condition, in these plots it artificially looks like the steady-state values
for the timecourses increase with concentration. This normalization is suspect,
however, because the maximum values reached are only around 2.5, whereas
previous datasets for the NBD-c126 mutant had shown the F/F0 values typically
reaching values between 3 and 4.  That, plus the fact that there is no good
reason for why the 20 nM cBid case would be more fluorescent other than that
the timecourses were cropped (e.g., there shouldn't be substantially
fluorescent buffer components for 20 nM cBid, for example) suggest that this
normalization is inappropriate.

Nevertheless, if we accept the normalization, we can fit the timecourses with
the single exponential function :math:`1 + F_{max} (1 - e^{-k_1 t})` and plot
the concentration dependence of the parameters:

.. ipython::

    In [8]: (k1, fmax, t0) = plot_normalized(bid_conc=20,
       ...: bid_conc_for_normalization=20, subtract_background=True,
       ...: do_fit=True, t0_val=None)

    @suppress
    In [15]: plt.savefig('_static/130911_c126_bax_titration_18.png')

.. image:: ../../_static/130911_c126_bax_titration_18.png
    :width: 6in

When we plot the fitted values for :math:`k_1`, it is clear that the insertion
process is slower at higher Bax:

.. ipython::

    In [16]: plt.figure()

    # We don't plot the fit values for 0 nM Bax
    In [17]: plt.plot(bax_concs[:-1], k1[:-1], linewidth=2, marker='o', color='r')

    @suppress
    In [18]: plt.xlabel('Bax (nM)')

    @suppress
    In [19]: plt.ylabel('$k_1$ value')

    @suppress
    In [20]: plt.savefig('_static/130911_c126_bax_titration_19.png')

.. image:: ../../_static/130911_c126_bax_titration_19.png
    :width: 6in

.. ipython::

    In [16]: plt.figure()

    # We don't plot the fit values for 0 nM Bax
    In [17]: plt.plot(bax_concs[:-1], fmax[:-1], linewidth=2, marker='o', color='r')

    @suppress
    In [18]: plt.xlabel('Bax (nM)')

    @suppress
    In [19]: plt.ylabel('$F_{max}$ value')

    @suppress
    In [20]: plt.savefig('_static/130911_c126_bax_titration_20.png')

.. image:: ../../_static/130911_c126_bax_titration_20.png
    :width: 6in

