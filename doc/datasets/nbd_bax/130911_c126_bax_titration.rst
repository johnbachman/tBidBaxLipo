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

Normalized by 0 cBid :math:`F_0`
--------------------------------

Now we look at the data normalized by the initial fluorescence values from the
0 cBid timecourses. Interestingly, though they are obviously noisy, the 0 nM
cBid timecourses themselves show that the F/F0 kinetics due to spontaneous Bax
insertion increases initially over the lower Bax concentrations, then seems to
plateau around 100 nM. There is a spontaneous increase in F/F0 in the 0 Bax
condition, though this likely represents a very small increase in absolute
terms that is exaggerated after normalization due to the fact that the F0 values
are small in this case. However, it might be worth subtracting this background
increase out from the raw fluorescence values before normalizing.

.. ipython::

    In [8]: plot_normalized_by_no_bid_f0(bid_conc=0)

    @suppress
    In [9]: plt.savefig('_static/130911_c126_bax_titration_4.png')

.. image:: ../../_static/130911_c126_bax_titration_4.png
    :width: 6in

With cBid added, from the shape of the curves it appears that there is a large
portion of the initial timecourses missing. That said, it appears that the
steady-state values for F/F0 increase for the lowest concentrations in the
titration, before plateauing around ~3.6 for 100-200 nM; however at this point
the kinetics start to slow down, even as the steady-state value appears to
remain constant.

.. ipython::

    In [10]: plot_normalized_by_no_bid_f0(bid_conc=20)

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_5.png')

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

    In [12]: (k1, fmax, t0) = plot_normalized_by_no_bid_f0_fits(bid_conc=20, \
       ....: t0_val=None)

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_5a.png')

.. image:: ../../_static/130911_c126_bax_titration_5a.png
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
irrelevant). Is this an artifact of the shape of the high concentration curves,
or does it reflect the possibility that the lower concentration wells were
pipetted first and the higher ones later?

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
    In [11]: plt.savefig('_static/130911_c126_bax_titration_5b.png')

.. image:: ../../_static/130911_c126_bax_titration_5b.png
    :width: 6in

On the other hand, the :math:`F_{max}` values show a hyperbolic increase,
plateauing around 200 nM:

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
    In [11]: plt.savefig('_static/130911_c126_bax_titration_5c.png')

.. image:: ../../_static/130911_c126_bax_titration_5c.png
    :width: 6in

Normalized by 0 nM cBid, fit with fixed :math:`t_0`
---------------------------------------------------

A potential concern about the results shown above is that the fitted values for
the parameters :math:`k_1` and :math:`F_{max}` may be influenced by the fitted
value for the parameter :math:`t_0`, which is allowed to vary for each curve.
Some of the fitted delays seem unrealistically long (e.g. ~800-900 sec for 12.5
and 25 nM Bax), and it also seems likely that a uniform delay for all wells due
to the software crash that Justin described would be large relative to the
smaller variable delays in pipetting the different wells.

To address this, we re-fit the data (normalized by the 0 nM cBid :math:`F_0`
values) while fixing :math:`t_0` to a constant value. We fix it to the average
over all the fitted values found in the previous step (ignoring the
fitted value for 0 Bax):

.. ipython::

    In [1]: mean_t0 = np.mean(t0[:-1]); print mean_t0

    In [2]: (k1, fmax, t0) = plot_normalized_by_no_bid_f0_fits(bid_conc=20, \
       ...: t0_val=mean_t0)

    @suppress
    In [11]: plt.savefig('_static/130911_c126_bax_titration_5d.png')

.. image:: ../../_static/130911_c126_bax_titration_5d.png
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
    In [11]: plt.savefig('_static/130911_c126_bax_titration_5e.png')

.. image:: ../../_static/130911_c126_bax_titration_5e.png
    :width: 6in

The :math:`F_{max}` values show a clear hyperbolic increase, plateauing around
200 nM:

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
    In [11]: plt.savefig('_static/130911_c126_bax_titration_5f.png')

.. image:: ../../_static/130911_c126_bax_titration_5f.png
    :width: 6in


Normalized by 20 nM cBid :math:`F_0`
------------------------------------

Alternatively, we can normalize the 20 nM cBid timecourses by their own initial value, which produces the following curves as a result:


.. ipython::

    In [12]: plot_normalized_by_20_bid_f0()

    @suppress
    In [13]: plt.savefig('_static/130911_c126_bax_titration_6.png')

.. image:: ../../_static/130911_c126_bax_titration_6.png
    :width: 6in

Here it looks like the steady-state values for the timecourses don't saturate
until fairly high concentrations of Bax, around 400 nM. One thing that makes
this normalization suspect, however, is that the maximum values reached are
only around 2.5, whereas previous datasets for the NBD-c126 mutant had shown
the F/F0 values typically reaching values between 3 and 4. That, plus the fact
that there is no good reason for why the 20 nM cBid case would be more
fluorescent other than that the timecourses were cropped (e.g., there shouldn't
be substantially fluorescent buffer components for 20 nM cBid, for example)
suggest that this normalization is inappropriate.

Nevertheless, if we accept the normalization, we can fit the timecourses with
the single exponential function :math:`1 + F_{max} (1 - e^{-k_1 t})` and plot
the concentration dependence of the parameters:

.. ipython::

    In [14]: (k1, fmax) = plot_normalized_by_20_bid_f0_fits()

    @suppress
    In [15]: plt.savefig('_static/130911_c126_bax_titration_7.png')

.. image:: ../../_static/130911_c126_bax_titration_7.png
    :width: 6in

The fits are fairly good, though not excellent; in particular they miss the
slight activation lag that is visible in the highest concentration curves. When we plot the fitted values for :math:`k_1`, it is clear that the insertion
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
    In [20]: plt.savefig('_static/130911_c126_bax_titration_8.png')

.. image:: ../../_static/130911_c126_bax_titration_8.png
    :width: 6in

For the :math:`F_{max}` values, we again see a hyperbolic increase with Bax
concentration. The principal difference with the results from the 0 nM cBid
normalization is that the :math:`F_{max}` values don't appear to saturate until
very high Bax concentrations (400-800 nM):

.. ipython::

    In [16]: plt.figure()

    # We don't plot the fit values for 0 nM Bax
    In [17]: plt.plot(bax_concs[:-1], fmax[:-1], linewidth=2, marker='o', color='r')

    @suppress
    In [18]: plt.xlabel('Bax (nM)')

    @suppress
    In [19]: plt.ylabel('$F_{max}$ value')

    @suppress
    In [20]: plt.savefig('_static/130911_c126_bax_titration_9.png')

.. image:: ../../_static/130911_c126_bax_titration_9.png
    :width: 6in

