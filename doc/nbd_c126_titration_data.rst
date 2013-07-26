Plots of Bax c126 NBD Titration Data
====================================

Raw data
--------

By Bax Concentration
~~~~~~~~~~~~~~~~~~~~

.. plot::
    :context:

    from matplotlib import pyplot as plt
    from matplotlib.font_manager import FontProperties
    import pickle

    data = pickle.load(open('../tbidbaxlipo/data/nbd_c126_titration_data.pck'))

    def do_legend():
        fontP = FontProperties()
        fontP.set_size=('small')
        ax = plt.gca()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        plt.legend(loc='upper left', prop=fontP, ncol=1, bbox_to_anchor=(1, 1),
                   fancybox=True, shadow=True)

    def plot_titration(bax_conc):
        timecourses = data.xs(bax_conc, axis=1, level='Bax')
        plt.figure()
        for lipo_conc in timecourses.columns:
            tc = timecourses[lipo_conc]
            plt.plot(tc.index.values, tc.values,
                     label='Lipo %.2f nM' % lipo_conc)
        plt.title('Liposome Titration for c126 Bax-NBD, %d nM' % bax_conc)
        do_legend()
        plt.show()

    for bax_conc in data.columns.levels[0]:
        plot_titration(bax_conc)


By Liposome Concentration
~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::
    :context:

    plt.close('all')

    def plot_titration(lipo_conc):
        timecourses = data.xs(lipo_conc, axis=1, level='Liposomes')
        plt.figure()
        for bax_conc in timecourses.columns:
            tc = timecourses[bax_conc]
            plt.plot(tc.index.values, tc.values, label='Bax %.2f nM' % bax_conc)
        plt.title('Bax c126 NBD Titration for %.2f nM Liposomes' % lipo_conc)
        do_legend()
        plt.show()

    for lipo_conc in data.columns.levels[1]:
        plot_titration(lipo_conc)

Background-Subtracted
---------------------

.. plot::
    :context:

    plt.close('all')

    def plot_bg_sub_titration(bax_conc):
        timecourses = data.xs(bax_conc, axis=1, level='Bax')
        plt.figure()
        bg = timecourses[0.]
        for lipo_conc in timecourses.columns:
            tc = timecourses[lipo_conc] - bg
            plt.plot(tc.index.values, tc.values,
                     label='Lipo %.2f nM' % lipo_conc)
        plt.title('Liposome Titration  for c126 Bax-NBD, %d nM' % bax_conc)
        do_legend()
        plt.show()

    for bax_conc in data.columns.levels[0]:
        plot_bg_sub_titration(bax_conc)

Single-Exponential Fits
-----------------------

.. plot::

    from matplotlib import pyplot as plt
    from matplotlib.font_manager import FontProperties
    import pickle
    from tbidbaxlipo.util import fitting
    import numpy as np

    data = pickle.load(open('../tbidbaxlipo/data/nbd_c126_titration_data.pck'))

    def plot_bg_sub_exp_fits(bax_conc):
        ks = []
        fmaxes = []
        timecourses = data.xs(bax_conc, axis=1, level='Bax')
        plt.figure()
        bg = timecourses[0.]
        time = np.array(bg.index.values, dtype='float')
        lipo_concs = np.array(timecourses.columns.values, dtype='float')
        for lipo_conc in lipo_concs:
            tc = timecourses[lipo_conc] - bg
            k = fitting.Parameter(np.log(2)/2300.)
            fmax = fitting.Parameter(3.)
            def single_exp(t):
                return fmax()*(1 - np.exp(-k()*t))
            fitting.fit(single_exp, [k, fmax], tc.values, time)
            ks.append(k())
            fmaxes.append(fmax())
            plt.plot(time, tc.values)
            plt.plot(time, np.array(map(single_exp, time)))
        plt.title('Single Exponential Fits to c126 Titration, Bax %d nM' % \
                  bax_conc)
        # Fit ks with powerlaw
        ks = np.array(ks)
        b = fitting.Parameter(1e-5)
        m = fitting.Parameter(1)
        def power_law(x): return b() * (x ** m())
        fitting.fit(power_law, [b, m], ks, lipo_concs)
        # Fit ks with hill func
        km = fitting.Parameter(1.)
        vmax = fitting.Parameter(1e-3)
        def hill_func(x): return (vmax() * x) / (km() + x)
        fitting.fit(hill_func, [km, vmax], ks, lipo_concs)
        # Plot data and fits
        plt.figure()
        plt.loglog(lipo_concs, ks, marker='o', color='b')
        plt.loglog(lipo_concs, map(power_law, lipo_concs), color='r',
                   label='Power')
        plt.loglog(lipo_concs, map(hill_func, lipo_concs), color='g',
                   label='Hill')
        plt.title('k vs. Liposome Concentration')
        print "b: %f" % b()
        print "m: %f" % m()
        plt.figure()
        plt.loglog(lipo_concs, np.array(fmaxes), marker='o')
        plt.title('Fmax vs. Liposome Concentration')
        plt.show()

    plt.ion()
    for bax_conc in data.columns.levels[0]:
        plot_bg_sub_exp_fits(bax_conc)
