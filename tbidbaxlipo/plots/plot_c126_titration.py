from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
import pickle
from tbidbaxlipo.util import fitting
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import Normalize
import pkgutil

# This allows me to get the data from any working directory
raw_data = pickle.loads(pkgutil.get_data('tbidbaxlipo.data',
                                     'nbd_c126_titration_data.pck'))

def plot_exponential_fits(data):
    plt.ion()
    bax_concs = np.array(data.columns.levels[0].values, dtype='float')
    lipo_concs = np.array(data.columns.levels[1].values[1:], dtype='float')

    k_matrix = np.zeros((len(bax_concs), len(lipo_concs)))
    fmax_matrix = np.zeros((len(bax_concs), len(lipo_concs)))

    for i, bax_conc in enumerate(bax_concs):
        plt.figure()
        #ks_for_lipo_concs = []
        #fmaxs_for_lipo_concs = []

        for j, lipo_conc in enumerate(lipo_concs):
            tc = data[(bax_conc, lipo_conc)]
            time = np.array(tc.index.values, dtype='float')
            k = fitting.Parameter(np.log(2)/2300.)
            fmax = fitting.Parameter(3.)
            def single_exp(t): return fmax()*(1 - np.exp(-k()*t))
            fitting.fit(single_exp, [k, fmax],
                        np.array(tc.values, dtype='float'), time)

            k_matrix[i, j] = k()
            fmax_matrix[i, j] = fmax()

            #ks_for_lipo_concs.append(k())
            #fmaxs_for_lipo_concs.append(fmax())

            plt.plot(time, tc.values)
            plt.plot(time, np.array(map(single_exp, time)))

        plt.title('Single Exponential Fits to c126 Titration, Bax %d nM' % \
                  bax_conc)

        ks_for_lipo_concs = k_matrix[i,:]
        fmaxs_for_lipo_concs = fmax_matrix[i,:]

        # Fit ks_for_lipo_concs with powerlaw
        b = fitting.Parameter(1e-5)
        m = fitting.Parameter(1)
        def power_law(x): return b() * (x ** m())
        fitting.fit(power_law, [b, m], ks_for_lipo_concs, lipo_concs)

        # Fit ks_for_lipo_concs with hill func
        km = fitting.Parameter(1.)
        vmax = fitting.Parameter(1e-3)
        def hill_func(x): return (vmax() * x) / (km() + x)
        fitting.fit(hill_func, [km, vmax], ks_for_lipo_concs, lipo_concs)

        # Plot data and fits
        plt.figure()
        plt.loglog(lipo_concs, ks_for_lipo_concs, marker='o', color='b')
        plt.loglog(lipo_concs, map(power_law, lipo_concs), color='r',
                   label='Power')
        plt.loglog(lipo_concs, map(hill_func, lipo_concs), color='g',
                   label='Hill')
        plt.title('k vs. Liposome Concentration')
        plt.figure()
        plt.loglog(lipo_concs, np.array(fmaxs_for_lipo_concs), marker='o')
        plt.title('Fmax vs. Liposome Concentration')
        plt.show()

        #k_matrix.append(ks_for_lipo_concs)
        #fmax_matrix.append(fmaxs_for_lipo_concs)

    return k_matrix

def background_subtract_data(data):
    """Returns a dataset with the no-liposome condition values subtracted."""
    bgsub_data = data.copy()
    bax_concs = data.columns.levels[0]
    lipo_concs = data.columns.levels[1]

    for bax_conc in data.columns.levels[0]:
        timecourses = data.xs(bax_conc, axis=1, level='Bax')
        bg = timecourses[0.]
        for lipo_conc in lipo_concs:
            bgsub_tc = timecourses[lipo_conc] - bg
            bgsub_data[(bax_conc, lipo_conc)] = bgsub_tc

    return bgsub_data

def do_legend():
    fontP = FontProperties()
    fontP.set_size=('small')
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.70, box.height])
    plt.legend(loc='upper left', prop=fontP, ncol=1, bbox_to_anchor=(1, 1),
               fancybox=True, shadow=True)

def plot_titration(data, outer_conc, outer_level):
    timecourses = data.xs(outer_conc, axis=1, level=outer_level)
    plt.figure()
    inner_level = timecourses.columns.name
    for inner_conc in timecourses.columns:
        tc = timecourses[inner_conc]
        plt.plot(tc.index.values, tc.values,
                 label='%s %.2f nM' % (inner_level, inner_conc))
    plt.title('Liposome/NBD-126-Bax titration, %s %.1f nM' %
              (outer_level, outer_conc))
    do_legend()
    plt.show()

# Functions used for making figures in docs
def plot_titrations_by_liposome_conc():
    for bax_conc in raw_data.columns.levels[0]:
        plot_titration(raw_data, bax_conc, 'Bax')

def plot_titrations_by_bax_conc():
    for lipo_conc in raw_data.columns.levels[1]:
        plot_titration(raw_data, lipo_conc, 'Liposomes')

def plot_insertion_rate_surface(k_matrix, data):
    lipo_concs = np.array(data.columns.levels[1].values[1:], dtype='float')
    bax_concs = np.array(data.columns.levels[0].values, dtype='float')
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(bax_concs, lipo_concs)
    ax = fig.gca()
    ax.plot_surface(X, Y, k_matrix.T, cmap=cm.coolwarm,
                    linewidth=1, rstride=1, cstride=1)
    ax.set_xlabel('[Bax] (nM)')
    ax.set_ylabel('[Lipos] (nM)')
    ax.set_zlabel('$k$ (sec$^{-1})$')
    plt.show()

