import numpy as np
from matplotlib import pyplot as plt
from preprocess_data import timecourse_wells, lipo_bg_wells, bgsub_wells, \
                            lipo_bg_layout, data_to_fit, bg_time
from tbidbaxlipo.util.plate_assay import plot_all, TIME, VALUE
from tbidbaxlipo.util import fitting

def plot_data():
    """Plots the data and various transformations of it."""
    # Timecourse wells
    plt.figure()
    plot_all(timecourse_wells)
    plt.title("Raw timecourses")

    # Lipo bg wells
    plt.figure()
    plot_all(lipo_bg_wells)
    plt.title("Lipo background wells")

    # Background-subtracted wells
    plt.figure()
    plot_all(bgsub_wells)
    plt.title("Background-subtracted wells")

    # Background-normalized wells
    plt.figure()
    for i, y in enumerate(data_to_fit):
        plt.plot(bg_time, y)
    plt.title('Background-normalized wells')

def plot_lipo_background(wells, layout):
    """Takes the lipo background well timecourses and plots the
    fluorescence values of the initial points as a function of lipo
    concentration. Shows that the liposome background fluorescence
    is strictly linear.
    """
    init_vals = np.zeros(len(layout.keys()))
    conc_list = np.zeros(len(layout.keys()))
    for i, cond_name in enumerate(layout):
        conc = float(cond_name.split(' ')[1])
        well_name = layout[cond_name][0]
        init_val = wells[well_name][VALUE][0]
        init_vals[i] = init_val
        conc_list[i] = conc

    # Fit the liposome background to a line
    print "Fitting liposome background"
    m = fitting.Parameter(1.)
    b = fitting.Parameter(1.)
    def linear(s):
        return m()*s + b()
    result = fitting.fit(linear, [m, b], init_vals, conc_list)

    plt.figure()
    plt.plot(conc_list, init_vals, marker='o', linestyle='', color='b')
    plt.plot(conc_list, linear(conc_list), color='r')
    plt.xlabel('Liposome concentration (mg/ml)')
    plt.ylabel('RFU')
    plt.title('Liposome background fluorescence at t=0')
    ax = plt.gca()
    ax.set_xscale('log')
    print "m: %f" % m()
    print "b: %f" % b()
    return [m(), b()]

if __name__ == '__main__':
    plt.ion()
    plot_data()
    plot_lipo_background(timecourse_wells, lipo_bg_layout)




