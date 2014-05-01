from tbidbaxlipo.util.plate_assay import *
import itertools
import pickle
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
from tbidbaxlipo.util import fitting
from matplotlib.font_manager import FontProperties

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict

layout = collections.OrderedDict([
        # Donor plus acceptor, competitor dilution
        ('cBid 2000 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A01', 'B01', 'C01']),
        ('cBid 1000 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A02', 'B02', 'C02']),
        ('cBid 500 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A03', 'B03', 'C03']),
        ('cBid 250 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A04', 'B04', 'C04']),
        ('cBid 125 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A05', 'B05', 'C05']),
        ('cBid 62.5 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A06', 'B06', 'C06']),
        ('cBid 31.3 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A07', 'B07', 'C07']),
        ('cBid 15.6 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A08', 'B08', 'C08']),
        ('cBid 7.8 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A09', 'B09', 'C09']),
        ('cBid 3.9 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A10', 'B10', 'C10']),
        ('cBid 1.95 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A11', 'B11', 'C11']),
        ('cBid 0 nM, cBid-488 20 nM, 1.55 nM Rho lipos', ['A12', 'B12', 'C12']),
        # No donor (to measure bleedthrough from acceptor and competitor)
        ('cBid 2000 nM, 1.55 nM Rho lipos', ['D01']),
        ('cBid 1000 nM, 1.55 nM Rho lipos', ['D02']),
        ('cBid 500 nM, 1.55 nM Rho lipos', ['D03']),
        ('cBid 250 nM, 1.55 nM Rho lipos', ['D04']),
        ('cBid 125 nM, 1.55 nM Rho lipos', ['D05']),
        ('cBid 62.5 nM, 1.55 nM Rho lipos', ['D06']),
        ('cBid 31.3 nM, 1.55 nM Rho lipos', ['D07']),
        ('cBid 15.6 nM, 1.55 nM Rho lipos', ['D08']),
        ('cBid 7.8 nM, 1.55 nM Rho lipos', ['D09']),
        ('cBid 3.9 nM, 1.55 nM Rho lipos', ['D10']),
        ('cBid 1.95 nM, 1.55 nM Rho lipos', ['D11']),
        ('cBid 0 nM, 1.55 nM Rho lipos', ['D12']),
        ('No donor', ['D01', 'D02', 'D03', 'D04', 'D05', 'D06', 'D07',
                           'D08', 'D09', 'D10', 'D11', 'D12']),
        # Donor with no acceptor, competitor dilution
        ('cBid 2000 nM, cBid-488 20 nM, 1.55 nM lipos', ['E01', 'F01', 'G01']),
        ('cBid 1000 nM, cBid-488 20 nM, 1.55 nM lipos', ['E02', 'F02', 'G02']),
        ('cBid 500 nM, cBid-488 20 nM, 1.55 nM lipos', ['E03', 'F03', 'G03']),
        ('cBid 250 nM, cBid-488 20 nM, 1.55 nM lipos', ['E04', 'F04', 'G04']),
        ('cBid 125 nM, cBid-488 20 nM, 1.55 nM lipos', ['E05', 'F05', 'G05']),
        ('cBid 62.5 nM, cBid-488 20 nM, 1.55 nM lipos', ['E06', 'F06', 'G06']),
        ('cBid 31.3 nM, cBid-488 20 nM, 1.55 nM lipos', ['E07', 'F07', 'G07']),
        ('cBid 15.6 nM, cBid-488 20 nM, 1.55 nM lipos', ['E08', 'F08', 'G08']),
        ('cBid 7.8 nM, cBid-488 20 nM, 1.55 nM lipos', ['E09', 'F09', 'G09']),
        ('cBid 3.9 nM, cBid-488 20 nM, 1.55 nM lipos', ['E10', 'F10', 'G10']),
        ('cBid 1.95 nM, cBid-488 20 nM, 1.55 nM lipos', ['E11', 'F11', 'G11']),
        ('cBid 0 nM, cBid-488 20 nM, 1.55 nM lipos', ['E12', 'F12', 'G12']),
        # No donor, no acceptor (to measure background from lipos)
        ('cBid 2000 nM, 1.55 nM lipos', ['H01']),
        ('cBid 1000 nM, 1.55 nM lipos', ['H02']),
        ('cBid 500 nM, 1.55 nM lipos', ['H03']),
        ('cBid 250 nM, 1.55 nM lipos', ['H04']),
        ('cBid 125 nM, 1.55 nM lipos', ['H05']),
        ('cBid 62.5 nM, 1.55 nM lipos', ['H06']),
        ('cBid 31.3 nM, 1.55 nM lipos', ['H07']),
        ('cBid 15.6 nM, 1.55 nM lipos', ['H08']),
        ('cBid 7.8 nM, 1.55 nM lipos', ['H09']),
        ('cBid 3.9 nM, 1.55 nM lipos', ['H10']),
        ('cBid 1.95 nM, 1.55 nM lipos', ['H11']),
        ('cBid 0 nM, 1.55 nM lipos', ['H12']),
        ('No donor no acceptor', ['H01', 'H02', 'H03', 'H04', 'H05', 'H06',
                                'H07', 'H08', 'H09', 'H10', 'H11', 'H12']),
    ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                 '140403_Bid_comp_timecourse.csv'))
preread_file = os.path.abspath(os.path.join(data_path,
                                 '140403_Bid_comp_preread.csv'))

ion()

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_wallac(timecourse_file)
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Preread wells
preread_wells = read_wallac(preread_file)

# Averages of raw timecourses across replicates
(timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)
"""Averages of raw timecourses."""

#
no_d_no_a_set = extract(
        ['cBid 2000 nM, 1.55 nM lipos',
        'cBid 1000 nM, 1.55 nM lipos',
        'cBid 500 nM, 1.55 nM lipos',
        'cBid 250 nM, 1.55 nM lipos',
        'cBid 125 nM, 1.55 nM lipos',
        'cBid 62.5 nM, 1.55 nM lipos',
        'cBid 31.3 nM, 1.55 nM lipos',
        'cBid 15.6 nM, 1.55 nM lipos',
        'cBid 7.8 nM, 1.55 nM lipos',
        'cBid 3.9 nM, 1.55 nM lipos',
        'cBid 1.95 nM, 1.55 nM lipos',
        'cBid 0 nM, 1.55 nM lipos'], timecourse_averages)
no_d_set = extract(
        ['cBid 2000 nM, 1.55 nM Rho lipos',
        'cBid 1000 nM, 1.55 nM Rho lipos',
        'cBid 500 nM, 1.55 nM Rho lipos',
        'cBid 250 nM, 1.55 nM Rho lipos',
        'cBid 125 nM, 1.55 nM Rho lipos',
        'cBid 62.5 nM, 1.55 nM Rho lipos',
        'cBid 31.3 nM, 1.55 nM Rho lipos',
        'cBid 15.6 nM, 1.55 nM Rho lipos',
        'cBid 7.8 nM, 1.55 nM Rho lipos',
        'cBid 3.9 nM, 1.55 nM Rho lipos',
        'cBid 1.95 nM, 1.55 nM Rho lipos',
        'cBid 0 nM, 1.55 nM Rho lipos'], timecourse_averages)

bg_all = extract(['No donor', 'No donor no acceptor'], timecourse_averages)
bg_all_sds = extract(['No donor', 'No donor no acceptor'], timecourse_stds)

d_a_set = extract(
        ['cBid 2000 nM, cBid-488 20 nM, 1.55 nM Rho lipos',
        'cBid 1000 nM, cBid-488 20 nM, 1.55 nM Rho lipos',
        'cBid 500 nM, cBid-488 20 nM, 1.55 nM Rho lipos',
        'cBid 250 nM, cBid-488 20 nM, 1.55 nM Rho lipos',
        'cBid 125 nM, cBid-488 20 nM, 1.55 nM Rho lipos',
        'cBid 62.5 nM, cBid-488 20 nM, 1.55 nM Rho lipos',
        'cBid 31.3 nM, cBid-488 20 nM, 1.55 nM Rho lipos',
        'cBid 15.6 nM, cBid-488 20 nM, 1.55 nM Rho lipos',
        'cBid 7.8 nM, cBid-488 20 nM, 1.55 nM Rho lipos',
        'cBid 3.9 nM, cBid-488 20 nM, 1.55 nM Rho lipos',
        'cBid 1.95 nM, cBid-488 20 nM, 1.55 nM Rho lipos',
        'cBid 0 nM, cBid-488 20 nM, 1.55 nM Rho lipos'], timecourse_averages)
d_set = extract(
        ['cBid 2000 nM, cBid-488 20 nM, 1.55 nM lipos',
        'cBid 1000 nM, cBid-488 20 nM, 1.55 nM lipos',
        'cBid 500 nM, cBid-488 20 nM, 1.55 nM lipos',
        'cBid 250 nM, cBid-488 20 nM, 1.55 nM lipos',
        'cBid 125 nM, cBid-488 20 nM, 1.55 nM lipos',
        'cBid 62.5 nM, cBid-488 20 nM, 1.55 nM lipos',
        'cBid 31.3 nM, cBid-488 20 nM, 1.55 nM lipos',
        'cBid 15.6 nM, cBid-488 20 nM, 1.55 nM lipos',
        'cBid 7.8 nM, cBid-488 20 nM, 1.55 nM lipos',
        'cBid 3.9 nM, cBid-488 20 nM, 1.55 nM lipos',
        'cBid 1.95 nM, cBid-488 20 nM, 1.55 nM lipos',
        'cBid 0 nM, cBid-488 20 nM, 1.55 nM lipos'], timecourse_averages)

"""
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
"""

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
    plot_all(d_set)
    title("Donor only")

    myfig()
    plot_all(d_a_set)
    title("Donor and acceptor")

    myfig()
    plot_all(no_d_no_a_set)
    title("No donor, no acceptor")

    myfig()
    plot_all(no_d_set)
    title("No donor")

    myfig()
    plot_all(bg_all, errors=bg_all_sds)
    title("Background averages")


    #myfig()
    #plot_all(bgsub_timecourses, errors=timecourse_stds, legend_position=0.7)
    #title("Averaged and bgsub")

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

def plot_fret(donor_acceptor, donor_only):
    # Plot each fret curve, F_DA/F_D
    figure(figsize=(11, 6))
    fret = collections.OrderedDict()
    for i, conc_str in enumerate(donor_acceptor.keys()):
        time = donor_acceptor[conc_str][TIME]
        f_d_plus_a = donor_acceptor[conc_str][VALUE]
        f_d = donor_only[donor_only.keys()[i]][VALUE]
        f = f_d_plus_a / f_d
        fret[conc_str] = []
        fret[conc_str].append(time)
        fret[conc_str].append(f)
        plot(time, f, label=conc_str)
        #f0 = fitting.Parameter(1)
        #k = fitting.Parameter(0.001)
        #fmin = fitting.Parameter(0.6)
        #def exp_decay(t):
        #    return f0()*np.exp(-k()*t) + fmin()
        #fitting.fit(exp_decay, [f0, k, fmin], f, time)
        #plot(time, exp_decay(time), marker='', color='gray')
    ylabel('$F_{D+A} / F_A$')
    xlabel('Time (sec)')
    title("FRET timecourses")
    fontP = FontProperties()
    fontP.set_size('small')
    ax = gca()
    box = ax.get_position()
    legend_position = 0.8
    ax.set_position([box.x0, box.y0, box.width * legend_position,
                     box.height])
    legend(loc='upper left', prop=fontP, ncol=1, bbox_to_anchor=(1, 1),
         fancybox=True, shadow=True)
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
    #plot_data()
    fret = plot_fret(d_a_set, d_set)
    #bid_only = plot_fit_donor_only()
    #fret = plot_fret_exp_fits(bid_only)
    #plot_fret_titration(fret, 0, 1)
    #plot_fret_titration(fret, -3, None)

