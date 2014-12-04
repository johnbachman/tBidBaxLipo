"""A script to parse the Bid-568 + DiD liposome FRET data directly from the
Excel spreadsheet sent by Justin Kale into Python data files."""

from openpyxl import load_workbook
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import collections
import itertools
from tbidbaxlipo.util import fitting

# Load the data
data_file = r'2014-04-29 - cBid-A568 FRET with DiD traced lipos - cBid unlab competition 1%w-w DiD.xlsx'
wb = load_workbook(data_file)

TIME_ROW_INDEX = 42
WELL_NAME_COL_INDEX = 0
FIRST_COL_INDEX = 1
FIRST_ROW_INDEX = 44
LAST_COL_INDEX = 154
LAST_ROW_INDEX = 140

# Load the first worksheet
sheet = wb.worksheets[0]

# Load the time coordinates
time = np.array([cell.value for cell in
                 sheet.rows[TIME_ROW_INDEX][FIRST_COL_INDEX:LAST_COL_INDEX]])

# FDA: WT cBid titration, DiD lipos, 10 nM Bid 568
fda = collections.OrderedDict([
         ('cBid 1000 nM', ['A1', 'B1', 'C1']),
         ('cBid 500 nM', ['A2', 'B2', 'C2']),
         ('cBid 250 nM', ['A3', 'B3', 'C3']),
         ('cBid 125 nM', ['A4', 'B4', 'C4']),
         ('cBid 62.5 nM', ['A5', 'B5', 'C5']),
         ('cBid 31.25 nM', ['A6', 'B6', 'C6']),
         ('cBid 15.6 nM', ['A7', 'B7', 'C7']),
         ('cBid 7.81 nM', ['A8', 'B8', 'C8']),
         ('cBid 3.91 nM', ['A9', 'B9', 'C9']),
         ('cBid 1.95 nM', ['A10', 'B10', 'C10']),
         ('cBid 0.98 nM', ['A11', 'B11', 'C11']),
         ('cBid 0 nM', ['A12', 'B12', 'C12']),
     ])

# FA: WT cBid titration, DiD lipos (no cBid-568)
fa = collections.OrderedDict([
         ('cBid 1000 nM', ['D1']),
         ('cBid 500 nM', ['D2']),
         ('cBid 250 nM', ['D3']),
         ('cBid 125 nM', ['D4']),
         ('cBid 62.5 nM', ['D5']),
         ('cBid 31.25 nM', ['D6']),
         ('cBid 15.6 nM', ['D7']),
         ('cBid 7.81 nM', ['D8']),
         ('cBid 3.91 nM', ['D9']),
         ('cBid 1.95 nM', ['D10']),
         ('cBid 0.98 nM', ['D11']),
         ('cBid 0 nM', ['D12'])
     ])
# FD: WT cBid titration, cBid-568, unlabeled lipos
fd = collections.OrderedDict([
         ('cBid 1000 nM', ['E1', 'F1', 'G1']),
         ('cBid 500 nM', ['E2', 'F2', 'G2']),
         ('cBid 250 nM', ['E3', 'F3', 'G3']),
         ('cBid 125 nM', ['E4', 'F4', 'G4']),
         ('cBid 62.5 nM', ['E5', 'F5', 'G5']),
         ('cBid 31.25 nM', ['E6', 'F6', 'G6']),
         ('cBid 15.6 nM', ['E7', 'F7', 'G7']),
         ('cBid 7.81 nM', ['E8', 'F8', 'G8']),
         ('cBid 3.91 nM', ['E9', 'F9', 'G9']),
         ('cBid 1.95 nM', ['E10', 'F10', 'G10']),
         ('cBid 0.98 nM', ['E11', 'F11', 'G11']),
         ('cBid 0 nM', ['E12', 'F12', 'G12']),
     ])
# BG: WT cBid titration, unlabeled lipos
bg = collections.OrderedDict([
         ('cBid 1000 nM', ['H1']),
         ('cBid 500 nM', ['H2']),
         ('cBid 250 nM', ['H3']),
         ('cBid 125 nM', ['H4']),
         ('cBid 62.5 nM', ['H5']),
         ('cBid 31.25 nM', ['H6']),
         ('cBid 15.6 nM', ['H7']),
         ('cBid 7.81 nM', ['H8']),
         ('cBid 3.91 nM', ['H9']),
         ('cBid 1.95 nM', ['H10']),
         ('cBid 0.98 nM', ['H11']),
         ('cBid 0 nM', ['H12'])
     ])

# The timecourses for each well will go in here
timecourse_wells = collections.OrderedDict()

# Iterate over the rows to get the individual trajectories
for row in sheet.rows[FIRST_ROW_INDEX:LAST_ROW_INDEX]:
    well_name = row[WELL_NAME_COL_INDEX].value
    timecourse = [cell.value for cell in row[FIRST_COL_INDEX:LAST_COL_INDEX]]
    timecourse_wells[well_name] = [time, timecourse]

### PLOTTING FUNCTIONS ### 

def plot_fd_final_values(wells, fd):
    (fd_avgs, fd_stds) = averages(wells, fd)
    plt.figure('FD')
    plot_all(fd_avgs)
    # Curiously, the overall level of the signal appears to be strongly
    # determined by the WT cBid concentration, but not in a monotonic way?
    bid_concs =  []
    final_val_avgs = []
    for conc_name in fd_avgs:
        bid_concs.append(float(conc_name.split(' ')[1]))
        final_val_avg = np.mean(fd_avgs[conc_name][VALUE][-10:])
        final_val_avgs.append(final_val_avg)
    plt.figure('FD values, final')
    plt.plot(bid_concs, final_val_avgs, marker='o')

def get_average_fa(wells, fa, plot=True):
    """Plot all FA (acceptor only, no donor) trajectories and then return the
    average across all of them. This average is used to do background subtraction
    for the FDA (donor + acceptor) trajectories.

    The fluorescence of the FA condition is not sensitive to the amount of
    unlabeled cBid present, so we can calculate our background vector by
    averaging across all of the conditions.

    Parameters
    ----------
    wells : collections.OrderedDict
        Dict containing the timecourses in each well.
    fa : collections.OrderedDict
        Dict mapping the concentration conditions to sets of wells for
        the FA condition.
    plot : boolean
        Whether to plot the timecourses.

    Return
    ------
    numpy.array containing the average FA timecourse.
    """

    # Calculate average background
    fa_well_names = itertools.chain(*fa.values())
    fa_wells = extract(fa_well_names, timecourse_wells)
    fa_avg = np.zeros((len(fa_wells.keys()), len(time)))
    for i, well_name in enumerate(fa_wells.keys()):
        fa_avg[i] = wells[well_name][VALUE]
    fa_avg = np.mean(fa_avg, axis=0)
    # Plot
    if plot:
        plt.figure('FA')
        plot_all(fa_wells)
        plt.plot(time, fa_avg, linewidth=3, color='k')
    # Return
    return fa_avg

def get_average_bg(wells, bg, plot=True):
    """Plot all BG (no donor, no acceptor) trajectories and then return the
    average across all of them. This average is used to do background subtraction
    for the FD (donor only, no acceptor) trajectories.

    The fluorescence of the BG condition is not sensitive to the amount of
    unlabeled cBid present, so we can calculate our background vector by
    averaging across all of the conditions.

    Parameters
    ----------
    wells : collections.OrderedDict
        Dict containing the timecourses in each well.
    bg : collections.OrderedDict
        Dict mapping the concentration conditions to sets of wells for
        the FA condition.
    plot : boolean
        Whether to plot the timecourses.

    Return
    ------
    numpy.array containing the average BG timecourse.
    """
    # Calculate average background
    bg_well_names = itertools.chain(*bg.values())
    bg_wells = extract(bg_well_names, timecourse_wells)
    bg_avg = np.zeros((len(bg_wells.keys()), len(time)))
    for i, well_name in enumerate(bg_wells.keys()):
        bg_avg[i] = wells[well_name][VALUE]
    bg_avg = np.mean(bg_avg, axis=0)
    # Plot
    if plot:
        plt.figure('BG')
        plot_all(bg_wells)
        plt.plot(time, bg_avg, linewidth=3, color='k')
    # Return
    return bg_avg

if __name__ == '__main__':
    from tbidbaxlipo.util.plate_assay import *
    plt.ion()

    plot_fd_final_values(timecourse_wells, fd)
    bg_avg = get_average_bg(timecourse_wells, bg, plot=True)
    fa_avg = get_average_fa(timecourse_wells, fa, plot=True)

    # Now, plot the curves, with fits

    for conc_name in fda.keys():
        plt.figure()
        # FDA
        for well_name in fda[conc_name]:
            fda_value = np.array(timecourse_wells[well_name][VALUE])
            # Subtract the FA background:
            fda_value = fda_value - fa_avg
            # Do the fit
            f0 = fitting.Parameter(3500.)
            fmax = fitting.Parameter(2000.)
            k = fitting.Parameter(1e-3)
            def fit_func(t):
                return fmax()*np.exp(-k()*t) + f0()
            res = fitting.fit(fit_func, [f0, fmax, k], fda_value, time)
            # Plot data and fit
            plt.plot(time, fda_value)
            plt.plot(time, fit_func(time), color='k', linewidth=2)
        # FD
        for well_name in fd[conc_name]:
            fd_value = np.array(timecourse_wells[well_name][VALUE])
            # Subtract the BG background:
            fd_value = fd_value - bg_avg
            # Do the fit
            f0 = fitting.Parameter(3500.)
            fmax = fitting.Parameter(2000.)
            k = fitting.Parameter(1e-3)
            def fit_func(t):
                return fmax()*np.exp(-k()*t) + f0()
            res = fitting.fit(fit_func, [f0, fmax, k], fd_value, time)
            # Plot data and fit
            plt.plot(time, fd_value)
            plt.plot(time, fit_func(time), color='k', linewidth=2)

# 1. Get residuals from all fits to estimate variance

# 2. Get f0, fmax, and k for all fits (fda and fd separately) to get estimates
#    of mean and variance of these parameters across replicates

# 3. Do hierarchical emcee fits for each concentration separately

# 4. Make predictive distributions over wells, and estimate equilibrium
#    FRET

# 5. Plot binding curve with associated error bars

