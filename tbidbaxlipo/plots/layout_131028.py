from tbidbaxlipo.util.plate_assay import plot_all, read_wallac, TIME, VALUE
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict

A_well_names = ['%s%.2d' % (row, col)
              for row in ['A', 'C', 'E','G','I','K','M','O']
              for col in (2 * np.arange(12)+1)]
B_well_names = ['%s%.2d' % (row, col)
              for row in ['B', 'D', 'F','H','J','L','N','P']
              for col in (2 * np.arange(12)+2)]

layout = collections.OrderedDict([
    ('Fluorescein 25 nM', A_well_names),
    ('Fluorescein 23.75 nM', B_well_names),
    ])

plt.ion()

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '131028_Fluorescein_384.csv'))

timecourse_wells = read_wallac(timecourse_file)
"""The raw (unnormalized) timecourses."""

A_timecourses = extract(A_well_names, timecourse_wells)
B_timecourses = extract(B_well_names, timecourse_wells)

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
    plt.hist(means, color='r')
    plt.title("Distribution of means for individual wells")
    plt.xlabel("Timecourse mean")
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

def plot_edge_effects(wells):
    # Now we plot across rows to look for edge effects
    wells_by_row = []
    for row in ['A', 'B', 'C','D','E','F','G','H','I','J','K','L','M','N','O','P']:
        row_wells = []
        for col in range(1,25):
            well_name = '%s%.2d' % (row, col)
            row_wells.append(well_name)
        wells_by_row.append(row_wells)

    wells_by_col = []
    for col in range(1, 25):
        col_wells = []
        for row in ['A', 'B', 'C','D','E','F','G','H','I','J','K','L','M','N','O','P']:
            well_name = '%s%.2d' % (row, col)
            col_wells.append(well_name)
        wells_by_col.append(col_wells)

    fig = plt.figure()
    for row in wells_by_row:
        row_vals = [np.mean(wells[well_name][VALUE])
                    for well_name in row]
        plt.plot(range(1,25), row_vals)
    plt.title('Edge effects vs. column')
    fontP = FontProperties()
    fontP.set_size('small')
    ax = fig.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.9, box.height])
    plt.legend(('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H','I','J','K','L','M','N',
                'O','P'), loc='upper left', prop=fontP, ncol=1,
                bbox_to_anchor=(1,1), fancybox=True, shadow=True)

    fig = plt.figure()
    for col in wells_by_col:
        col_vals = [np.mean(wells[well_name][VALUE])
                    for well_name in col]
        plt.plot(col_vals)
    plt.title('Edge effects vs. row')
    fontP = FontProperties()
    fontP.set_size('small')
    ax = fig.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    plt.legend([str(num) for num in range(1,25)],
                loc='upper left', prop=fontP, ncol=2,
                bbox_to_anchor=(1,1), fancybox=True, shadow=True)

