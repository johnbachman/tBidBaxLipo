Plots of Bax c126 NBD Titration Data
====================================

Raw data
--------

By Bax Concentration
~~~~~~~~~~~~~~~~~~~~

.. plot::

    from matplotlib import pyplot as plt
    from matplotlib.font_manager import FontProperties
    import pickle

    data = pickle.load(open('../tbidbaxlipo/data/nbd_c126_titration_data.pck'))

    def plot_titration(bax_conc):
        timecourses = data.xs(bax_conc, axis=1, level='Bax')
        plt.figure()
        for lipo_conc in timecourses.columns:
            tc = timecourses[lipo_conc]
            plt.plot(tc.index.values, tc.values,
                     label='Lipo %.2f nM' % lipo_conc)
        plt.title('Liposome Titration for c126 Bax-NBD, %d nM' % bax_conc)

        fontP = FontProperties()
        fontP.set_size=('small')
        ax = plt.gca()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        plt.legend(loc='upper left', prop=fontP, ncol=1, bbox_to_anchor=(1, 1),
                   fancybox=True, shadow=True)
        plt.show()

    for bax_conc in data.columns.levels[0]:
        plot_titration(bax_conc)

By Liposome Concentration
~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::

    from matplotlib import pyplot as plt
    from matplotlib.font_manager import FontProperties
    import pickle

    data = pickle.load(open('../tbidbaxlipo/data/nbd_c126_titration_data.pck'))

    def plot_titration(lipo_conc):
        timecourses = data.xs(lipo_conc, axis=1, level='Liposomes')
        plt.figure()
        for bax_conc in timecourses.columns:
            tc = timecourses[bax_conc]
            plt.plot(tc.index.values, tc.values, label='Bax %.2f nM' % bax_conc)
        plt.title('Bax c126 NBD Titration for %.2f nM Liposomes' % lipo_conc)

        fontP = FontProperties()
        fontP.set_size=('small')
        ax = plt.gca()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        plt.legend(loc='upper left', prop=fontP, ncol=1, bbox_to_anchor=(1, 1),
                   fancybox=True, shadow=True)
        plt.show()

    for lipo_conc in data.columns.levels[1]:
        plot_titration(lipo_conc)

