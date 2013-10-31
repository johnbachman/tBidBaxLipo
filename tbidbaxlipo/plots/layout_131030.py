from tbidbaxlipo.util.plate_assay import *
import itertools
import pickle
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
from tbidbaxlipo.util import fitting

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict

layout = collections.OrderedDict([
        ('Fluorescein 100 nM', ['A1', 'B1', 'C1', 'D1', 'E1']),
        ('Fluorescein 50 nM', ['A2', 'B2', 'C2', 'D2', 'E2']),
        ('Fluorescein 25 nM', ['A3', 'B3', 'C3', 'D3', 'E3']),
        ('Fluorescein 12.5 nM', ['A4', 'B4', 'C4', 'D4', 'E4']),
        ('Fluorescein 6.25 nM', ['A5', 'B5', 'C5', 'D5', 'E5']),
        ('Fluorescein 3.125 nM', ['A6', 'B6', 'C6', 'D6', 'E6']),
        ('Fluorescein 1.56 nM', ['A7', 'B7', 'C7', 'D7', 'E7']),
        ('Fluorescein 0.781 nM', ['A8', 'B8', 'C8', 'D8', 'E8']),
        ('Fluorescein 0.391 nM', ['A9', 'B9', 'C9', 'D9', 'E9']),
        ('Fluorescein 0.195 nM', ['A10', 'B10', 'C10', 'D10', 'E10']),
        ('Fluorescein 0.098 nM', ['D11', 'B11', 'C11', 'D11', 'E11']),
        ('Fluorescein 0 nM', ['D12', 'B12', 'C12', 'D12', 'E12']),
        ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)

num_concs = 12
num_reps = 5
num_timepoints = 20

ion()

def get_wells(filename):
    timecourse_file = os.path.abspath(os.path.join(data_path, filename))
                                     #'131030_Fluorescein_Flex_Med_6read.txt'))
                                     #'131022_Fluorescein_01sec.csv'))
                                     #'131021_Fluorescein_03sec.csv'))

    # Timecourse wells
    timecourse_wells = read_flex(timecourse_file)

    # Averages of raw timecourses across replicates
    (timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)
    """Averages of raw timecourses."""

    # Fill in data matrix
    data_matrix = np.zeros((num_concs, num_reps, num_timepoints))
    conc_list = []
    for conc_index, conc_str in enumerate(layout.keys()):
        conc_wells = layout[conc_str]
        conc_list.append(float(conc_str.split(' ')[1]))
        for rep_index, well_name in enumerate(conc_wells):
            data_matrix[conc_index, rep_index, :] = \
                    timecourse_wells[well_name][VALUE]
    conc_list = np.array(conc_list)

    # Get mean across all timepoints
    means = np.mean(data_matrix, axis=2)
    stds = np.std(data_matrix, axis=2)

    # Take logs of data
    log_concs = np.log10(conc_list)
    log_data = np.log10(data_matrix)
    log_means = np.mean(log_data, axis=2)
    log_stds = np.std(log_data, axis=2)
    return (timecourse_wells, timecourse_averages, timecourse_stds,
            conc_list, means, stds, log_concs, log_means, log_stds)

def plot_data(timecourse_wells, timecourse_averages, timecourse_stds):
    def myfig():
        figure(figsize=(11, 6))

    myfig()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    myfig()
    plot_all(timecourse_averages, errors=timecourse_stds)
    title("Raw timecourses, averaged across replicates")

def plot_cvs(means, stds, log_concs):
    # Plot CV at each concentration for the replicates
    figure()
    cvs = (stds / means) * 100.
    for rep_index in range(cvs.shape[1]):
        plot(log_concs, cvs[:,rep_index], marker='o',
             label='Row %d' % (rep_index + 1))
    title('CV of %s reads vs. [Fluorescein]' % num_timepoints)
    xlabel('log10([Fluorescein]) (nM)')
    ylabel('CV (%)')
    legend(loc='upper right')

def plot_dilution_series(conc_list, means, stds, log_concs, log_means, log_stds):
    # Plot each dilution series separately
    figure()
    for rep_index in range(means.shape[1]):
        errorbar(conc_list, means[:,rep_index], yerr=stds[:,rep_index],
                 label='Row %d' % (rep_index + 1))
    legend(loc='lower right')
    xlabel('[Fluorescein] (nM)')
    ylabel('RFU')
    title('Fluorescein dilution series')

    # And now on a log scale
    figure()
    for rep_index in range(log_means.shape[1]):
        errorbar(log_concs, log_means[:,rep_index], yerr=log_stds[:, rep_index],
                 label='Row %d' % (rep_index + 1))
    xlabel('log10([Fluorescein]) (nM)')
    ylabel('log10(RFU)')
    legend(loc='lower right')
    title('Fluorescein dilution series, log-log plot')

def plot_bar_plot(means, stds):
    # Bar plot of replicates showing error bars
    figure()
    width = 1 / float((means.shape[1]+1))
    for rep_index in range(means.shape[1]):
        bar(np.arange(num_concs)+(width*rep_index), means[:,rep_index], width,
                 yerr=stds[:,rep_index])

def plot_fits(means, log_means, conc_list, log_concs):
    # Ignore any NaNs
    # Take mean of all dilution series
    dilution_means = np.mean(means, axis=1)
    nans = np.isnan(dilution_means)
    dilution_means = dilution_means[nans==False]
    dilution_stds = np.std(means[nans==False], axis=1)

    # Take mean of all dilution series
    log_dilution_means = np.mean(log_means[nans==False], axis=1)
    log_dilution_stds = np.std(log_means[nans==False], axis=1)

    conc_list_no_nan = conc_list[nans == False]
    log_concs_no_nan = log_concs[nans == False]


    # Fit with line
    m = fitting.Parameter(900)
    b = fitting.Parameter(0)
    def linear(x):
        return m()*x + b()
    fitting.fit(linear, [m, b], dilution_means, conc_list_no_nan)
    print m()
    print b()

    figure()
    errorbar(conc_list_no_nan, dilution_means, yerr=dilution_stds,
            color='r', linestyle='')
    plot(conc_list_no_nan, linear(conc_list_no_nan), color='r', label='Linear')
    #plot(conc_list_no_nan, quenching_quad(conc_list_no_nan), color='g',
    #     label='Dimerization')
    #plot(conc_list, andersson(conc_list), color='b', label='Andersson et al.')
    #plot(conc_list, power_law(conc_list), color='m', label='Power law')
    xlabel('[Fluorescein] (nM)')
    ylabel('RFU')
    title('Fits to Fluorescein titration')
    legend(loc='lower right')

    figure()
    errorbar(log_concs_no_nan, log_dilution_means, yerr=log_dilution_stds,
             color='r', linestyle='')
    plot(log_concs_no_nan, np.log10(linear(conc_list_no_nan)),
            color='r', label='Linear')
    legend(loc='lower right')
    xlabel('log10([Fluorescein]) (nM)')
    ylabel('log10(RFU)')
    title('Fits to Fluorescein titration, log-log')

    return fitting.r_squared(dilution_means, linear(conc_list_no_nan))


if __name__ == '__main__':
    for filename in ['131030_Fluorescein_Flex_High_22read.txt',
                     '131030_Fluorescein_Flex_Med_6read.txt',
                     '131030_Fluorescein_Flex_Med_100read.txt']:
        [timecourse_wells, timecourse_averages, timecourse_stds,
         conc_list, means, stds, log_concs, log_means, log_stds] = \
                        get_wells(filename)

        plot_data(timecourse_wells, timecourse_averages, timecourse_stds)
        plot_bar_plot(means, stds)
        plot_cvs(means, stds, log_concs)
        plot_dilution_series(conc_list, means, stds, log_concs, log_means,
                             log_stds)
        r_squared = plot_fits(means, log_means, conc_list, log_concs)
        print "R^2: %f" % r_squared
