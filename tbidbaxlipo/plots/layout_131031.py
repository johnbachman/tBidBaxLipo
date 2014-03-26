from tbidbaxlipo.util.plate_assay import *
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
import numpy as np

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict

fluorescein_well_names = ['%s%d' % (row, col)
              for row in ['A', 'B', 'C','D','E','F','G']
              for col in range(1, 13)]
buffer_well_names = ['%s%d' % (row, col)
              for row in ['H']
              for col in range(1, 13)]

plt.ion()

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                            #'131030_Fluorescein_Flex_50uL_83nM_med_12read.txt'))
                            #'131031_BotRead_Fluor25nM_50uL.txt'))
                            '131031_TopRead_Fluor25nM_50uL.txt'))

timecourse_wells = read_flexstation_flex(timecourse_file)
"""The raw (unnormalized) timecourses."""

def get_means_sds_cvs(wells):
    # Calculate the mean, SD, and CV for the timecourse in each well
    means = []
    sds = []
    cvs = []
    for well in wells.keys():
        value_arr = np.array(wells[well][VALUE])
        #if np.isnan(value_arr).any():
         #   import pdb; pdb.set_trace()
        value_arr = value_arr[np.isnan(value_arr) == False]
        means.append(np.mean(value_arr))
        sds.append(np.std(value_arr))
    for i, mean in enumerate(means):
        cvs.append(100 * (sds[i] / float(mean)))
    return (means, sds, cvs)

def plot(wells):
    (means, sds, cvs) = get_means_sds_cvs(wells)

    plt.figure()
    plot_all(wells)
    plt.title("Timecourses, raw")

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

if __name__ == '__main__':
    plot_all(timecourse_wells)
    timecourse_wells = truncate_timecourses(timecourse_wells, 12)
    fluorescein_wells = extract(fluorescein_well_names, timecourse_wells)
    plot_all(timecourse_wells)
    [means, sds, cvs] = get_means_sds_cvs(fluorescein_wells)
    plot(fluorescein_wells)
    print_statistics(fluorescein_wells)


