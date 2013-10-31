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
        ('cBid 50 nM + 15.5 nM liposomes', ['E01', 'F01', 'G01']),
        ('cBid 50 nM + 7.75 nM liposomes', ['E02', 'F02', 'G02']),
        ('cBid 50 nM + 3.875 nM liposomes', ['E03', 'F03', 'G03']),
        ('cBid 50 nM + 1.938 nM liposomes', ['E04', 'F04', 'G04']),
        ('cBid 50 nM + 0.969 nM liposomes', ['E05', 'F05', 'G05']),
        ('cBid 50 nM + 0.484 nM liposomes', ['E06', 'F06', 'G06']),
        ('cBid 50 nM + 0.242 nM liposomes', ['E07', 'F07', 'G07']),
        ('cBid 50 nM + 0.121 nM liposomes', ['E08', 'F08', 'G08']),
        ('cBid 50 nM + 0.061 nM liposomes', ['E09', 'F09', 'G09']),
        ('cBid 50 nM + 0.030 nM liposomes', ['E10', 'F10', 'G10']),
        ('cBid 50 nM + 0.015 nM liposomes', ['E11', 'F11', 'G11']),
        ('cBid 50 nM + 0 nM liposomes', ['E12', 'F12', 'G12']),
        ('15.5 nM liposomes', ['H01']),
        ('7.75 nM liposomes', ['H02']),
        ('3.875 nM liposomes', ['H03']),
        ('1.938 nM liposomes', ['H04']),
        ('0.969 nM liposomes', ['H05']),
        ('0.484 nM liposomes', ['H06']),
        ('0.242 nM liposomes', ['H07']),
        ('0.121 nM liposomes', ['H08']),
        ('0.061 nM liposomes', ['H09']),
        ('0.030 nM liposomes', ['H10']),
        ('0.015 nM liposomes', ['H11']),
        ('0 nM liposomes', ['H12']),
    ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                 '131016_Bid488_RhodPE_fret.csv'))

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

signal_set = extract(
       ['cBid 50 nM + 15.5 nM liposomes',
        'cBid 50 nM + 7.75 nM liposomes',
        'cBid 50 nM + 3.875 nM liposomes',
        'cBid 50 nM + 1.938 nM liposomes',
        'cBid 50 nM + 0.969 nM liposomes',
        'cBid 50 nM + 0.484 nM liposomes',
        'cBid 50 nM + 0.242 nM liposomes',
        'cBid 50 nM + 0.121 nM liposomes',
        'cBid 50 nM + 0.061 nM liposomes',
        'cBid 50 nM + 0.030 nM liposomes',
        'cBid 50 nM + 0.015 nM liposomes',
        'cBid 50 nM + 0 nM liposomes'], timecourse_averages)

background_set = extract(
       ['15.5 nM liposomes',
        '7.75 nM liposomes',
        '3.875 nM liposomes',
        '1.938 nM liposomes',
        '0.969 nM liposomes',
        '0.484 nM liposomes',
        '0.242 nM liposomes',
        '0.121 nM liposomes',
        '0.061 nM liposomes',
        '0.030 nM liposomes',
        '0.015 nM liposomes',
        '0 nM liposomes'], timecourse_averages)

bgsub_timecourses = subtract_background_set(signal_set, background_set)

# Initial timepoints
#initial_wells_tc = get_first_points_by_well(timecourse_wells)

# Normalized timecourses
#norm_wells = get_normalized_well_timecourses(
#        timecourse_wells, initial_wells_tc, final_well_avgs)

# Normalized timecourses, averaged
#(norm_averages, norm_stds) = averages(norm_wells, layout)

# First timepoint shifted to 0 (better for fitting)
#reset_norm_subset = reset_first_timepoint_to_zero(norm_averages)

def plot_data():
    def myfig():
        figure(figsize=(11, 6))

    myfig()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    myfig()
    plot_all(timecourse_averages, errors=timecourse_stds, legend_position=0.7)
    title("Raw timecourses, averaged")

    myfig()
    plot_all(bgsub_timecourses, errors=timecourse_stds, legend_position=0.7)
    title("Averaged and bgsub")

def plot_fit_donor_only():
    # Fit the background (donor-only) to a line
    m = fitting.Parameter(0)
    b = fitting.Parameter(26000)
    def linear(x):
        return m()*x + b()
    t = timecourse_averages['cBid 50 nM + 0 nM liposomes'][TIME]
    x = timecourse_averages['cBid 50 nM + 0 nM liposomes'][VALUE]
    fitting.fit(linear, [m, b], x, t)
    figure()
    plot(t, x)
    plot(t, linear(t))
    xlabel('Time (sec)')
    ylabel('RFU')
    title('Linear fit of Bid-only control')
    return linear(t)

def plot_fret_exp_fits(donor_only):
    # Plot each fret curve, F_DA/F_D
    figure()
    fret = collections.OrderedDict()
    f_d = bgsub_timecourses['cBid 50 nM + 0 nM liposomes'][VALUE]
    for i, conc_str in enumerate(bgsub_timecourses.keys()):
        time = bgsub_timecourses[conc_str][TIME]
        f_d_plus_a = bgsub_timecourses[conc_str][VALUE]
        f = f_d_plus_a / donor_only
        fret[conc_str] = []
        fret[conc_str].append(time)
        fret[conc_str].append(f)
        plot(time, f, marker='o', linestyle='')

        f0 = fitting.Parameter(1)
        k = fitting.Parameter(0.001)
        fmin = fitting.Parameter(0.6)
        def exp_decay(t):
            return f0()*np.exp(-k()*t) + fmin()
        fitting.fit(exp_decay, [f0, k, fmin], f, time)
        plot(time, exp_decay(time), marker='', color='gray')
    ylabel('$F_{D+A} / F_A$')
    xlabel('Time (sec)')
    title("FRET timecourses")
    return fret

def plot_fret_titration(fret, timepoint_start, timepoint_end):
    # Calculate fret at each concentration
    lipo_concs = []
    fret_0 = []
    for i, conc_str in enumerate(fret.keys()):
        #if i == 0:
        #    continue
        lipo_concs.append(float(conc_str.split(' ')[4]))
        f_0 = np.mean(fret[conc_str][VALUE][timepoint_start:timepoint_end])
        fret_0.append(f_0)
    lipo_concs = np.array(lipo_concs)
    fret_0 = np.array(fret_0)

    # Plot curve as up from 0
    figure()
    plot(lipo_concs, 1 - fret_0, marker='o', color='r', label='Data')
    ylabel('$F_{D+A} / F_D$')
    xlabel('[Liposomes] (nM)')
    title('FRET titration')
    # Fit with binding curve
    kd = fitting.Parameter(3.)
    fmax = fitting.Parameter(0.25)
    def binding_eqn(x):
        return (fmax() * x) / (kd() + x)
    fitting.fit(binding_eqn, [kd, fmax], 1 - fret_0, lipo_concs)
    plot(lipo_concs, binding_eqn(lipo_concs), color='g', label='Hill')
    # Fit with exp curve
    kd = fitting.Parameter(0.6/3.)
    fmax = fitting.Parameter(0.60)
    def exp_eqn(x):
        return fmax() * (1 - np.exp(-kd()*x))
    fitting.fit(exp_eqn, [kd, fmax], 1 - fret_0, lipo_concs)
    plot(lipo_concs, binding_eqn(lipo_concs), color='m', label='Exponential')
    # Fit with NSB curve
    nsb_kd = fitting.Parameter(3.)
    nsb_fmax = fitting.Parameter(0.25)
    nsb_nsb = fitting.Parameter(0.001)
    def nsb_eqn(x):
        return ((nsb_fmax() * x) / (nsb_kd() + x)) + nsb_nsb()*x
    fitting.fit(nsb_eqn, [nsb_kd, nsb_fmax, nsb_nsb], 1 - fret_0, lipo_concs)
    plot(lipo_concs, nsb_eqn(lipo_concs), color='b', label='Nonsat. binding')
    legend(loc='lower right')

    # Fit with quadratic curve
    lipid_concs = lipo_concs * 83775.
    fmax = fitting.Parameter(0.25)
    kd = fitting.Parameter(3)
    def quad(atot):
        btot = 50.
        return fmax() * (
               (1 / (2 * btot)) *
               (atot + btot + kd() -
                np.sqrt((atot + btot + kd())**2 - (4 * atot * btot)))
               )
    fitting.fit(quad, [kd, fmax], 1 - fret_0, lipid_concs)
    figure()
    plot(lipid_concs, 1 - fret_0, color='r', marker='o')
    ylabel('$F_{D+A} / F_D$')
    xlabel('[Lipid] (nM)')
    title('FRET titration')
    plot(lipid_concs, quad(lipid_concs), color='b', label='Quadratic')
    legend(loc='lower right')


if __name__ == '__main__':
    plot_data()
    bid_only = plot_fit_donor_only()
    fret = plot_fret_exp_fits(bid_only)
    plot_fret_titration(fret, 0, 1)
    plot_fret_titration(fret, -3, None)

