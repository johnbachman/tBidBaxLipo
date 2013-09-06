from tbidbaxlipo.util.plate_assay import plot_all, read_wallac, TIME, VALUE
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
import numpy as np

all_well_names = ['%s%.2d' % (row, col)
              for row in ['A', 'B', 'C','D','E','F','G','H']
              for col in range(1, 13)]

plt.ion()

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '130906_wallac_dispenser_50ul.csv'))

timecourse_wells = read_wallac(timecourse_file)
"""The raw (unnormalized) timecourses."""

def get_means_sds_cvs(wells):
    # Calculate the mean, SD, and CV for the timecourse in each well
    means = []
    sds = []
    cvs = []
    for well in wells.keys():
        means.append(np.mean(wells[well][VALUE]))
        sds.append(np.std(wells[well][VALUE]))
    for i, mean in enumerate(means):
        cvs.append(100 * (sds[i] / float(mean)))
    return (means, sds, cvs)

def plot(wells):
    (means, sds, cvs) = get_means_sds_cvs(wells)

    plt.figure()
    plot_all(wells)
    plt.title("20uL timecourses, raw")

    plt.figure()
    plt.hist(means, color='r')
    plt.title("Distribution of means for individual wells")
    plt.xlabel("Timecourse mean")
    plt.ylabel("Count")

    plt.figure()
    plt.hist(sds, color='r')
    plt.title("Distribution of SDs for individual wells")
    plt.xlabel("Timecourse standard deviation")
    plt.ylabel("Count")

    plt.figure()
    plt.hist(cvs, color='r')
    plt.title("Distribution of CVs for individual wells")
    plt.xlabel("Timecourse coeff. of variation")
    plt.ylabel("Count")

def print_statistics(wells):
    (means, sds, cvs) = get_means_sds_cvs(wells)

    # Now, find the pipetting error
    global_mean = np.mean(means)
    global_sd = np.std(means)
    global_cv = 100 * (global_sd / global_mean)

    print "Mean across dispensed wells: %f" % global_mean
    print "SD across dispensed wells:   %f" % global_sd
    print "CV across dispensed wells:   %f" % global_cv

def plot_edge_effects():
    # Now we plot across rows to look for edge effects
    wells_by_row = []
    for row in ['A', 'B', 'C','D','E','F','G','H']:
        row_wells = []
        for col in range(1,13):
            well_name = '%s%.2d' % (row, col)
            row_wells.append(well_name)
        wells_by_row.append(row_wells)

    plt.figure()
    for row in wells_by_row:
        row_vals = [np.mean(timecourse_wells[well_name][VALUE])
                    for well_name in row]
        plt.plot(range(1,13), row_vals)
        plt.legend(('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'), loc='lower right')

