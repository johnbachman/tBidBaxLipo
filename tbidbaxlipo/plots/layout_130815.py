import numpy as np
from tbidbaxlipo.util.plate_assay import *
import tbidbaxlipo.data
import os

layout = collections.OrderedDict([
    ('Bax 4000 nM', ['H01', 'H03']),
    ('Bax 0 nM', ['H02', 'H04'])
    ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                               '130815_Baxa5_Timecourses.csv'))
triton_file = os.path.abspath(os.path.join(data_path,
                                           '130815_Baxa5_Triton.csv'))

timecourse_wells = read_wallac(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = reset_well_times(timecourse_wells, offset=1.)

# Initial timepoints
initial_wells_tc = get_first_points_by_well(timecourse_wells)

# Post-Triton values
final_wells = read_wallac(triton_file)
final_well_avgs = get_repeat_averages_by_well(final_wells)

# Normalized curves
norm_wells = get_normalized_well_timecourses(
        timecourse_wells, 8200., final_well_avgs)
"""Timecourses normalized to min/max (initial/Triton) values."""

(norm_averages, norm_stds) = averages(norm_wells, layout)
"""Timecourses normalized and then averaged."""

data = to_dataframe(norm_averages, norm_stds)
"""pandas DataFrame of data normalized, averaged data."""

def main():
    """Plots the data and various transformations of it."""
    ion()

    # Timecourse wells
    figure()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    # Normalized wells
    figure()
    plot_all(norm_wells)
    title("Normalized timecourses")

    # Normalized and averaged
    figure()
    plot_all(norm_averages)
    title("Normalized timecourses, averaged")

