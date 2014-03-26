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


"""
layout = collections.OrderedDict([
        ('Fluorescein 500 nM', ['C01', 'D01', 'E01']),
        ('Fluorescein 250 nM', ['C02', 'D02', 'E02']),
        ('Fluorescein 125 nM', ['C03', 'D03', 'E03']),
        ('Fluorescein 62.5 nM', ['C04', 'D04', 'E04']),
        ('Fluorescein 31.25 nM', ['C05', 'D05', 'E05']),
        ('Fluorescein 15.63 nM', ['C06', 'D06', 'E06']),
        ('Fluorescein 7.81 nM', ['C07', 'D07', 'E07']),
        ('Fluorescein 3.91 nM', ['C08', 'D08', 'E08']),
        ('Fluorescein 1.95 nM', ['C09', 'D09', 'E09']),
        ('Fluorescein 0.98 nM', ['C10', 'D10', 'E10']),
        ('Fluorescein 0.49 nM', ['C11', 'D11', 'E11']),
        ('Fluorescein 0.25 nM', ['C12', 'D12', 'E12']),
        #('Fluorescein 0 nM', ['G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07',
        #                      'G08', 'G09']),
        ])

layout = collections.OrderedDict([
        ('Fluorescein 1000 nM', ['A01', 'B01', 'C01']),
        ('Fluorescein 500 nM', ['A02', 'B02', 'C02']),
        ('Fluorescein 250 nM', ['A03', 'B03', 'C03']),
        ('Fluorescein 125 nM', ['A04', 'B04', 'C04']),
        ('Fluorescein 62.5 nM', ['A05', 'B05', 'C05']),
        ('Fluorescein 31.25 nM', ['A06', 'B06', 'C06']),
        ('Fluorescein 15.63 nM', ['A07', 'B07', 'C07']),
        ('Fluorescein 7.81 nM', ['A08', 'B08', 'C08']),
        ('Fluorescein 3.91 nM', ['A09', 'B09', 'C09']),
        ('Fluorescein 1.95 nM', ['A10', 'B10', 'C10']),
        ('Fluorescein 0.98 nM', ['A11', 'B11', 'C11']),
        ('Fluorescein 0.49 nM', ['A12', 'B12', 'C12']),
        ])
"""

layout = collections.OrderedDict([
        ('Fluorescein 200 nM', ['D01', 'E01', 'F01', 'G01', 'H01']),
        ('Fluorescein 100 nM', ['D02', 'E02', 'F02', 'G02', 'H02']),
        ('Fluorescein 50 nM', ['D03', 'E03', 'F03', 'G03', 'H03']),
        ('Fluorescein 25 nM', ['D04', 'E04', 'F04', 'G04', 'H04']),
        ('Fluorescein 12.5 nM', ['D05', 'E05', 'F05', 'G05', 'H05']),
        ('Fluorescein 6.25 nM', ['D06', 'E06', 'F06', 'G06', 'H06']),
        ('Fluorescein 3.125 nM', ['D07', 'E07', 'F07', 'G07', 'H07']),
        ('Fluorescein 1.56 nM', ['D08', 'E08', 'F08', 'G08', 'H08']),
        ('Fluorescein 0.781 nM', ['D09', 'E09', 'F09', 'G09', 'H09']),
        ('Fluorescein 0.391 nM', ['D10', 'E10', 'F10', 'G10', 'H10']),
        ('Fluorescein 0.195 nM', ['D11', 'E11', 'F11', 'G11', 'H11']),
        #('Fluorescein 0.098 nM', ['D12', 'E12', 'F12', 'G12', 'H12']),
        ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                 '131023_Fluorescein_1sec.csv'))
                                 #'131022_Fluorescein_01sec.csv'))
                                 #'131021_Fluorescein_03sec.csv'))

ion()

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_wallac(timecourse_file)
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Averages of raw timecourses across replicates
(timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)
"""Averages of raw timecourses."""

num_concs = 11
num_reps = 5
num_timepoints = 30
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


def plot_data():
    def myfig():
        figure(figsize=(11, 6))

    myfig()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    myfig()
    plot_all(timecourse_averages, errors=timecourse_stds)
    title("Raw timecourses, averaged across replicates")

def plot_cvs():
    # Plot CV at each concentration for the three replicates
    figure()
    cvs = (stds / means) * 100.
    for rep_index in range(cvs.shape[1]):
        plot(log_concs, cvs[:,rep_index], marker='o',
             label='Row %d' % (rep_index + 1))
    #plot(log_concs, cvs[:,1], marker='o', color='g', label='Row D')
    #plot(log_concs, cvs[:,2], marker='o', color='b', label='Row E')
    title('CV of %s reads vs. [Fluorescein]' % num_timepoints)
    xlabel('log10([Fluorescein]) (nM)')
    ylabel('CV (%)')
    legend(loc='upper right')

def plot_dilution_series():
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
        errorbar(log_concs, log_means[:,rep_index], yerr=log_stds[:,rep_index],
                 label='Row %d' % (rep_index + 1))
    xlabel('log10([Fluorescein]) (nM)')
    ylabel('log10(RFU)')
    legend(loc='lower right')
    title('Fluorescein dilution series, log-log plot')

def plot_bar_plot():
    # Bar plot of replicates showing error bars
    figure()
    width = 1 / float((means.shape[1]+1))
    for rep_index in range(means.shape[1]):
        bar(np.arange(num_concs)+(width*rep_index), means[:,rep_index], width,
                 yerr=stds[:,rep_index])

def plot_fits():
    # Take mean of all dilution series
    dilution_means = np.mean(means, axis=1)
    dilution_stds = np.std(means, axis=1)

    # Fit with line
    m = fitting.Parameter(900)
    b = fitting.Parameter(0)
    def linear(x):
        return m()*x + b()
    fitting.fit(linear, [m], dilution_means, conc_list)
    print b()
    # Fit with dimerization func
    kf = fitting.Parameter(1e-3)
    kr = fitting.Parameter(1e2)
    fmax = fitting.Parameter(900)
    def quenching_quad(x):
        return fmax()*((-kr() + np.sqrt(kr()**2 + 4*kf()*kr()*x)) /
                       (2 * kf())) + 35.
    fitting.fit(quenching_quad, [kf, kr, fmax], dilution_means, conc_list)

    # Fit with Andersson et al. func
    beta = fitting.Parameter(1000)
    ext_coeff = fitting.Parameter(70000*1e-9)
    k = fitting.Parameter(1)
    def andersson(x):
        #return beta()*x * np.exp(-k()*x) + 40
        return beta() * (1 - 10**(-ext_coeff()*x)) * np.exp(-k()*x)
    fitting.fit(andersson, [beta, ext_coeff, k], dilution_means, conc_list)

    # Fit with power law
    power_b = fitting.Parameter(2.)
    power_m = fitting.Parameter(1.)
    def power_law(x):
        return power_b()*(x ** power_m())
    fitting.fit(power_law, [power_b, power_m], dilution_means, conc_list)

    # Fit with power law
    log_b = fitting.Parameter(2.)
    log_m = fitting.Parameter(1.)
    def log_line(x):
        return log_m()*x + log_b()
    fitting.fit(log_line, [log_m, log_b], np.log10(dilution_means), log_concs)

    figure()
    errorbar(conc_list, dilution_means, yerr=dilution_stds, color='r',
             linestyle='')
    plot(conc_list, linear(conc_list), color='r', label='Linear')
    plot(conc_list, quenching_quad(conc_list), color='g', label='Dimerization')
    plot(conc_list, andersson(conc_list), color='b', label='Andersson et al.')
    plot(conc_list, power_law(conc_list), color='m', label='Power law')
    xlabel('[Fluorescein] (nM)')
    ylabel('RFU')
    title('Fits to Fluorescein titration')
    legend(loc='lower right')

    figure()
    errorbar(log_concs, log_means[:, 0], yerr=log_stds[:,0], color='k',
             linestyle='')
    plot(log_concs, np.log10(linear(conc_list)),
            color='r', label='Linear')
    plot(log_concs, np.log10(quenching_quad(conc_list)),
            color='g', label='Dimerization')
    plot(log_concs, np.log10(andersson(conc_list)),
            color='b', label='Andersson et al.')
    plot(log_concs, np.log10(power_law(conc_list)),
            color='m', label='Power law')
    plot(log_concs, log_line(log_concs),
            color='k', label='Log line')

    legend(loc='lower right')
    xlabel('log10([Fluorescein]) (nM)')
    ylabel('log10(RFU)')
    title('Fits to Fluorescein titration, log-log')

if __name__ == '__main__':
    plot_bar_plot()
    plot_cvs()
    plot_dilution_series()
    plot_fits()
    sys.exit()

    plot_data()

    figure()
    dilution_index = range(12)
    m = fitting.Parameter(1500)
    r = fitting.Parameter(0.9)
    def compound_error(n):
        return m()*r()**n
    dilution_ratio0 = means[:,0] / conc_list
    dilution_ratio1 = means[:,1] / conc_list
    dilution_ratio2 = means[:,2] / conc_list
    plot(dilution_index, dilution_ratio0, marker='o', linestyle='', color='r')
    plot(dilution_index, dilution_ratio1, marker='o', linestyle='', color='g')
    plot(dilution_index, dilution_ratio2, marker='o', linestyle='', color='b')

    fitting.fit(compound_error, [m, r], dilution_ratio0, dilution_index)
    plot(dilution_index, compound_error(dilution_index), color='r')
