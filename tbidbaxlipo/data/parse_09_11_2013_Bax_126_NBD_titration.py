"""
A script to parse the NBD fluorescence data from the Excel spreadsheet sent by
Justin Kale into Python data files.
"""

from openpyxl import load_workbook
import numpy as np
import pandas as pd
from itertools import product
import pickle

from matplotlib import pyplot as plt
from tbidbaxlipo.util import fitting
from texttable import Texttable
from matplotlib.font_manager import FontProperties
import tbidbaxlipo.data
import os
import sys

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
data_file = os.path.abspath(os.path.join(data_path,
        '09-11-2013 - Bax 126C NBD titration 10 nM lipos, 20 nM cbid.xlsx'))

fontP = FontProperties()
fontP.set_size('small')

# Zero-indexed
FIRST_COL_INDEX = 4
LAST_COL_INDEX = 204 # 336
FIRST_ROW_INDEX = 45
LAST_ROW_INDEX = 362
NUM_CONCS = 8

wb = load_workbook(data_file)

# Load the first worksheet
sheet0 = wb.worksheets[0]

# Load the time coordinates
time = np.array([cell.value for cell in
                 sheet0.rows[43][FIRST_COL_INDEX:LAST_COL_INDEX]])

# Load the liposome concentrations
bax_concs = np.array([cell.value for cell in
                      sheet0.columns[1][45:53]])


def get_data_from_titration(row_start, num_replicates):
    # Get data column by column as a list of lists
    data = []
    for conc_index in range(NUM_CONCS):
        for rep_index in range(num_replicates):
            row_arr = []
            data_row = sheet0.rows[row_start + conc_index +
                    (rep_index * NUM_CONCS)][FIRST_COL_INDEX:LAST_COL_INDEX]
            # Iterate and fill the array of column values
            for cell in data_row:
                if (cell.value == None):
                    row_arr.append(np.nan)
                else:
                    row_arr.append(cell.value)
            data.append(row_arr)
    return np.array(data)

# Get the data from each sheet
cbid0nm = get_data_from_titration(FIRST_ROW_INDEX, 1)
cbid20nm = get_data_from_titration(53, 3)
data_matrix = np.concatenate((cbid0nm, cbid20nm), axis=0)

# Build the column index
col_tuples = list(product([0.], bax_concs, [1]))
col_tuples += list(product([20.], bax_concs, [1, 2, 3]))

col_index = pd.MultiIndex.from_tuples(col_tuples,
                                      names=('Bid', 'Bax', 'Replicate'))

# Build the dataset
data = pd.DataFrame(data_matrix.T, index=time, columns=col_index)

norm_data_matrix = [tc / float(tc[0]) for tc in data_matrix]
norm_data_matrix = np.array(norm_data_matrix)
norm_data = pd.DataFrame(norm_data_matrix.T, index=time, columns=col_index)

def r_squared(y, yhat):
    ybar = np.mean(y)
    ssres = np.sum((y - yhat) ** 2)
    #ssreg = np.sum((yhat - ybar) ** 2)
    sstot = np.sum((y - ybar) ** 2)
    return 1 - (ssres / sstot)

def plot_f0():
    plt.figure()

    # Plot f0 values for no cBid
    f0_no_bid = []
    for conc_index in range(NUM_CONCS):
        conc_data = data[0][bax_concs[conc_index]][1].values
        f0_no_bid.append(conc_data[0])
    # (fit nonzero f0 values to a line)
    m = fitting.Parameter(100)
    b = fitting.Parameter(500)
    def linear(concs):
        return m() * concs + b()
    fitting.fit(linear, [m, b], f0_no_bid[:-1], bax_concs[:-1])
    plt.plot(bax_concs[:], f0_no_bid[:], marker='o',
                 color='g', label='0 cBid, $F_0$ data', linewidth=2)
    plt.plot(bax_concs, linear(bax_concs), label='0 cBid, fit', color='m',
             linewidth=2)

    # Print table of results
    tt = Texttable()
    tt.header(['[Bax]', 'F0'])
    tt.add_rows(reversed(zip(bax_concs, f0_no_bid)), header=False)
    print "========= F0 values, 0 cBid ========="
    print
    print tt.draw()
    print
    print "F0 linear fit, no Bid"
    print "m: %s" % m()
    print "b: %s" % b()
    print "R^2: %s" % r_squared(f0_no_bid, linear(bax_concs))
    print

    # Plot f0 values for 20 nM cBid
    f0_means = []
    f0_stds = []
    for conc_index in range(NUM_CONCS):
        conc_data = data[20][bax_concs[conc_index]][:].values
        conc_mean = np.mean(conc_data, axis=1)
        conc_std = np.std(conc_data, axis=1)
        f0_means.append(conc_mean[0])
        f0_stds.append(conc_std[0])
    # (fit nonzero f0 values to a line)
    f0_cvs = np.array(f0_stds) / np.array(f0_means)
    m = fitting.Parameter(100)
    b = fitting.Parameter(500)
    def linear(concs):
        return m() * concs + b()
    fitting.fit(linear, [m, b], f0_means[:-1], bax_concs[:-1])
    plt.errorbar(bax_concs[:], f0_means[:], yerr=f0_stds[:] / np.sqrt(3),
                 linewidth=2, color='b', label='20 nM cBid, $F_0$ data')
    plt.plot(bax_concs, linear(bax_concs), color='r', linewidth=2,
             label='20 nM cBid, fit')

    # Print table of results
    tt = Texttable()
    tt.header(['[Bax]', 'F0 (mean)', 'CV (%)'])
    tt.add_rows(reversed(zip(bax_concs, f0_means, f0_cvs)), header=False)
    print "========= F0 values, 20 nM cBid ========="
    print
    print tt.draw()
    print
    print "Fitting F0 values to y = mx + b..."
    print "NOTE: Ignoring F0 for Bax = 0 in fit"
    print "m: %s" % m()
    print "b: %s" % b()
    print "R^2: %s" % r_squared(f0_means[:-1], linear(bax_concs[:-1]))

    plt.legend(loc='lower right')
    plt.title('$F_0$ values')
    plt.show()

def plot_raw_timecourses(bid_conc=0):
    # 20 nM cBid
    plt.figure()
    for conc_index in range(NUM_CONCS):
        conc_data = data[bid_conc][bax_concs[conc_index]][:].values
        time = data[bid_conc][bax_concs[conc_index]][:].index
        plt.errorbar(time, np.mean(conc_data, axis=1),
                 yerr=np.std(conc_data, axis=1) / np.sqrt(3),
                 label='Bax %s nM' % bax_concs[conc_index])
    box = plt.gca().get_position()
    plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plt.legend(loc='upper left', prop=fontP, ncol=1, bbox_to_anchor=(1, 1),
           fancybox=True, shadow=True)
    plt.title('Raw data, %s nM cBid' % bid_conc)
    plt.show()

def plot_normalized(bid_conc=0, bid_conc_for_normalization=0,
                    subtract_background=True, do_fit=False,
                    t0_val=None):
    plt.figure()
    k1_arr = []
    fmax_arr = []
    t0_arr = []
    for conc_index in range(NUM_CONCS):
        time = np.array(data[bid_conc][bax_concs[conc_index]][:].index.values,
                        dtype='float')
        conc_data = data[bid_conc][bax_concs[conc_index]][:].values
        f0 = data[bid_conc_for_normalization] \
                 [bax_concs[conc_index]][1].values[0]
        if subtract_background: # from 0 Bax timecourse
            bg = data[bid_conc][0][:].mean(axis=1).values
            conc_data = np.subtract(conc_data.T, bg).T
            f0 = f0 - bg[0]
        conc_data = conc_data / float(f0)
        # Plots the fits if desired
        if do_fit:
            fmax = fitting.Parameter(3.)
            k1 = fitting.Parameter(0.0005)
            def one_exp(t):
                return 1. + fmax() * (1 - np.exp(-k1() * (t + t0())))
            # Check if we're supposed to fit the t0 parameter or not
            if t0_val is None:
                t0 = fitting.Parameter(556)
                fitting.fit(one_exp, [fmax, k1, t0],
                            np.mean(conc_data, axis=1), time)
            else:
                t0 = fitting.Parameter(t0_val)
                fitting.fit(one_exp, [fmax, k1], np.mean(conc_data, axis=1),
                            time)

            k1_arr.append(k1())
            fmax_arr.append(fmax())
            t0_arr.append(t0())
            # Plot
            plt.plot(time, one_exp(time), color='k')
            plt.plot(time, np.mean(conc_data, axis=1),
                     label='Bax %s nM' % bax_concs[conc_index])
        else:
            plt.errorbar(time, np.mean(conc_data, axis=1),
                         yerr=np.std(conc_data, axis=1) / np.sqrt(3),
                         label='Bax %s nM' % bax_concs[conc_index])

    box = plt.gca().get_position()
    plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plt.legend(loc='upper left', prop=fontP, ncol=1,
               bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
    plt.title(r'Normalized ($F/F_0$) by %s nM cBid $F_0$' % \
              bid_conc_for_normalization)
    plt.show()

    if do_fit:
        return (k1_arr, fmax_arr, t0_arr)
    else:
        return None

if __name__ == '__main__':
    pass
