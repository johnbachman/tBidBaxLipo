from tbidbaxlipo.util.plate_assay import plot_all, read_wallac, TIME, VALUE
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
import string

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict

all_well_names = ['%s%.2d' % (row, col)
              for row in ['A', 'B', 'C','D','E','F','G','H','I','J','K','L',
                          'M','N','O','P']
              for col in range(1, 24)]

plt.ion()

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
#timecourse_file = os.path.abspath(os.path.join(data_path,
#                                        '140113_Fluorescein_384.csv'))
#timecourse_file = os.path.abspath(os.path.join(data_path,
#                                        '140113_Fluorescein_384_flipped.csv'))
#timecourse_file = os.path.abspath(os.path.join(data_path,
#                                        '140204_Fluorescein_384.csv'))
#timecourse_file = os.path.abspath(os.path.join(data_path,
#                                        '140219_ANTS_plt3881.csv'))
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '140221_ANTS_plt3686_edge_flipped.csv'))
#timecourse_file = os.path.abspath(os.path.join(data_path,
#                                        '140204_Fluorescein_384_flipped.csv'))

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

def get_mean_plate_map(wells, num_rows=8, num_cols=12):
    row_names = [string.uppercase[i] for i in range(num_rows)]
    col_names = ['%.2d' % num for num in range(1, num_cols+1)]
    plate_map = np.zeros((num_rows, num_cols))

    for i, row in enumerate(row_names):
        for j, col in enumerate(col_names):
            well_name = row + col
            plate_map[i, j] = np.mean(wells[well_name][VALUE])

    return plate_map

def plot(wells):
    (means, sds, cvs) = get_means_sds_cvs(wells)

    plt.figure()
    plt.hist(means, color='r')
    plt.title("Distribution of means for individual wells")
    plt.xlabel("Timecourse mean")
    plt.ylabel("Count")

    plt.figure()
    plt.hist(cvs, color='r')
    plt.title("Distribution of CVs for individual wells")
    plt.xlabel("Timecourse CV")
    plt.ylabel("Count")

    plate_map = get_mean_plate_map(wells)
    plt.figure()
    plt.pcolor(plate_map)
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.show()

def print_statistics(wells):
    (means, sds, cvs) = get_means_sds_cvs(wells)

    # Now, find the pipetting error
    global_mean = np.mean(means)
    global_sd = np.std(means)
    global_cv = 100 * (global_sd / global_mean)

    print "Mean across dispensed wells: %f" % global_mean
    print "SD across dispensed wells:   %f" % global_sd
    print "CV across dispensed wells:   %f" % global_cv

def plot_edge_effects(wells, num_rows=8, num_cols=12):
    # Now we plot across rows to look for edge effects
    wells_by_row = []
    row_names = [string.uppercase[i] for i in range(num_rows)]
    for row in row_names:
        row_wells = []
        for col in range(1,num_cols+1):
            well_name = '%s%.2d' % (row, col)
            row_wells.append(well_name)
        wells_by_row.append(row_wells)

    wells_by_col = []
    for col in range(1, num_cols+1):
        col_wells = []
        for row in row_names:
            well_name = '%s%.2d' % (row, col)
            col_wells.append(well_name)
        wells_by_col.append(col_wells)

    fig = plt.figure()
    for row in wells_by_row:
        row_vals = [np.mean(wells[well_name][VALUE])
                    for well_name in row]
        plt.plot(range(1,num_cols+1), row_vals)
    plt.title('Edge effects vs. column')
    fontP = FontProperties()
    fontP.set_size('small')
    ax = fig.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.9, box.height])
    plt.legend(row_names, loc='upper left', prop=fontP, ncol=1,
                bbox_to_anchor=(1,1), fancybox=True, shadow=True)
    plt.xlabel('Column index')
    plt.ylabel('RFU')

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
    plt.legend([str(num) for num in range(1,num_cols+1)],
                loc='upper left', prop=fontP, ncol=2,
                bbox_to_anchor=(1,1), fancybox=True, shadow=True)
    plt.xlabel('Row index')
    plt.ylabel('RFU')

