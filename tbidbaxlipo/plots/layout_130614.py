from tbidbaxlipo.util.plate_assay import *
import itertools
import pickle
import collections
import sys
import os
import tbidbaxlipo.data

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict

layout = collections.OrderedDict([
        #('Bax 0 nM',  ['A12', 'B12', 'C12']),
        ('Bax 11 nM', ['A11', 'B11', 'C11']),
        ('Bax 16 nM', ['A10', 'B10', 'C10']),
        ('Bax 25 nM', ['A09', 'B09', 'C09']),
        ('Bax 37 nM',   ['A08', 'B08', 'C08']),
        ('Bax 55 nM',  ['A07', 'B07', 'C07']),
        ('Bax 83 nM',  ['A06', 'B06', 'C06']),
        ('Bax 124 nM',  ['A05', 'B05', 'C05']),
        ('Bax 187 nM',  ['A04', 'B04', 'C04']),
        ('Bax 280 nM', ['A03', 'B03', 'C03']),
        ('Bax 420 nM', ['A02', 'B02', 'C02']),
        ('Bax 630 nM', ['A01', 'B01', 'C01']),
        ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '130614_Bax_43C_Timecourse.csv'))
triton_file = os.path.abspath(os.path.join(data_path,
                                        '130614_Bax_43C_Triton.csv'))

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_wallac(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Initial timepoints
initial_wells_tc = get_first_points_by_well(timecourse_wells)

# Post-Triton values
final_wells = read_wallac(triton_file)
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

# First timepoint shifted to 0 (better for fitting)
reset_norm_means = reset_first_timepoint_to_zero(norm_averages)
"""Timecourses normalized, averaged, then with first point shifted to t = 0."""
reset_norm_sds = reset_first_timepoint_to_zero(norm_stds)

# Pandas dataframe
df = to_dataframe(reset_norm_means, reset_norm_sds)
"""Pandas DataFrame version of reset_norm_means/sds"""

# Pore timecourses
pores = get_average_pore_timecourses(reset_norm_means)
"""Average pores derived by taking -log(1 - data)."""

def main():
    """Plots the data and various transformations of it."""
    ion()

    # Timecourse wells
    figure()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    # Averages of raw timecourses across replicates
    figure()
    plot_all(timecourse_averages, errors=timecourse_stds)
    title("Raw timecourses, averaged")

    # Normalized timecourses
    figure()
    plot_all(norm_wells)
    title("Normalized timecourses")

    # Normalized timecourses, averaged
    figure()
    plot_all(norm_averages, errors=norm_stds)
    title("Normalized timecourses, averaged")

    # First timepoint shifted to 0 (better for fitting)
    figure()
    plot_all(reset_norm_means)
    title("Norm., Avg., Reset to t = 0")

    # Pore timecourses
    figure()
    plot_all(pores)
    title("Avg. pores per liposome")
