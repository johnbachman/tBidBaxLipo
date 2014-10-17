from tbidbaxlipo.util.plate_assay import *
import itertools
import pickle
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import fitting
from tbidbaxlipo.plots import titration_fits
import pymc
from scipy import stats

# All wells have 100 nM cBid and 100 nM Bax
layout = collections.OrderedDict([
        ('Lipos 21.7 nM, Pre-inc',  ['C1', 'D1']),
        ('Lipos 10.9 nM, Pre-inc',  ['C2', 'D2']),
        ('Lipos 5.4 nM, Pre-inc',  ['C3', 'D3']),
        ('Lipos 2.7 nM, Pre-inc',  ['C4', 'D4']),
        ('Lipos 1.4 nM, Pre-inc',  ['C5', 'D5']),
        ('Lipos 0.68 nM, Pre-inc',  ['C6', 'D6']),
        ('Lipos 0.34 nM, Pre-inc',  ['C7', 'D7']),
        ('Lipos 0.17 nM, Pre-inc',  ['C8', 'D8']),
        ('Lipos 0.085 nM, Pre-inc',  ['C9', 'D9']),
        ('Lipos 0.042 nM, Pre-inc',  ['C10', 'D10']),
        ('Lipos 0.021 nM, Pre-inc',  ['C11', 'D11']),
        ('Lipos 0 nM, Pre-inc',  ['C12', 'D12']),
        ('Lipos 21.7 nM',  ['E1', 'F1']),
        ('Lipos 10.9 nM',  ['E2', 'F2']),
        ('Lipos 5.4 nM',  ['E3', 'F3']),
        ('Lipos 2.7 nM',  ['E4', 'F4']),
        ('Lipos 1.4 nM',  ['E5', 'F5']),
        ('Lipos 0.68 nM',  ['E6', 'F6']),
        ('Lipos 0.34 nM',  ['E7', 'F7']),
        ('Lipos 0.17 nM',  ['E8', 'F8']),
        ('Lipos 0.085 nM',  ['E9', 'F9']),
        ('Lipos 0.042 nM',  ['E10', 'F10']),
        ('Lipos 0.021 nM',  ['E11', 'F11']),
        ('Lipos 0 nM',  ['E12', 'F12']),
        ('ANTS Lipos 12.1 nM, Bid/Bax', ['G1', 'G2']),
        ('ANTS Lipos 12.1 nM', ['G3', 'G4']),
        ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
preinc_file = os.path.abspath(os.path.join(data_path,
                                        '141016_Bax_depletion_preincubation.txt'))
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '141016_Bax_depletion_timecourse.txt'))
triton_file = os.path.abspath(os.path.join(data_path,
                                        '141016_Bax_depletion_triton.txt'))

# Define a few sets of wells
release_ctrls = ['ANTS Lipos 12.1 nM, Bid/Bax',
                 'ANTS Lipos 12.1 nM']
release_ctrl_wells = ['G1', 'G2', 'G3', 'G4']

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Preincubation wells
preinc_wells = read_flexstation_kinetics(preinc_file)

# Timecourse wells
timecourse_wells = read_flexstation_kinetics(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Initial timepoints
initial_wells_tc = get_first_points_by_well(timecourse_wells)

# Post-Triton values
final_wells = read_flexstation_kinetics(triton_file)
final_wells = extract(wells_to_read, final_wells)
final_well_avgs = get_repeat_averages_by_well(final_wells)

# Averages of raw timecourses across replicates
(timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)
"""Averages of raw timecourses."""

# Normalized timecourses
norm_wells = get_normalized_well_timecourses(
        timecourse_wells, initial_wells_tc, final_well_avgs)
"""Timecourses normalized to min/max (initial/Triton) values."""

# Normalized timecourses, averaged
(norm_averages, norm_stds) = averages(norm_wells, layout)
"""Timecourses normalized and then averaged."""

# Get background average
#background = norm_averages['Bax 0 nM'][VALUE]

# Normalized and background subtracted
#bgsub_norm_wells = subtract_background(norm_wells, background)

# Normalized, background subtracted, averaged
#(bgsub_norm_averages, bgsub_norm_stds) = averages(bgsub_norm_wells, layout)

# First timepoint shifted to 0 (better for fitting)
#reset_bgsub_means = reset_first_timepoint_to_zero(bgsub_norm_averages)
"""Timecourses normalized, BG-subtracted, averaged, then with first point
shifted to t = 0."""
#reset_bgsub_sds = reset_first_timepoint_to_zero(bgsub_norm_stds)

# Pandas dataframe
#df = to_dataframe(bgsub_norm_averages, bgsub_norm_stds)
"""Pandas DataFrame version of reset_norm_means/sds"""

# Pore timecourses
#pores = get_average_pore_timecourses(bgsub_norm_averages)
"""Average pores derived by taking -log(1 - data)."""

def plot_data():
    """Plots the data and various transformations of it."""
    ion()
    # ANTS release ctrls
    ctrls = extract(release_ctrl_wells, preinc_wells)
    plot_all(ctrls)
    return
    # Timecourse wells
    figure()
    plot_all(timecourse_wells)
    title("Raw timecourses")
    return
    # Averages of raw timecourses across replicates
    figure()
    plot_all(timecourse_averages, errors=timecourse_stds)
    title("Raw timecourses, averaged")

    # Normalized timecourses
    figure()
    plot_all(norm_wells)
    title("Normalized timecourses")

    # Normalized timecourses, background-subtracted
    figure()
    plot_all(bgsub_norm_wells)
    title("Normalized, BG-subtracted timecourses")

    # Normalized timecourses, averaged
    figure()
    plot_all(norm_averages, errors=norm_stds)
    title("Normalized timecourses, averaged")

    # Normalized timecourses, background subtracted, averaged
    figure()
    plot_all(bgsub_norm_averages, errors=norm_stds)
    title("Normalized timecourses, BG-subtracted, averaged")

    # First timepoint shifted to 0 (better for fitting)
    figure()
    plot_all(reset_bgsub_means)
    title("Norm., BG-sub, avg., Reset to t = 0")

    # Pore timecourses
    figure()
    plot_all(pores)
    title("Avg. pores per liposome")

if __name__ == '__main__':
    # First, plot the ANTS controls to show that reaction should be at eq
    plot_data()
    # Second, big plot, maybe 3 x 4 subplots, showing 
