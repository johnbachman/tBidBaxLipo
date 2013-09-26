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
        ('Bax 0 nM, tBid 0 nM', ['E01']),
        ('Bax 0 nM, tBid 0.24 nM', ['E02']),
        ('Bax 0 nM, tBid 2.35 nM', ['E03']),
        ('Bax 0 nM, tBid 23.5 nM', ['E04']),
        #('Bax 0 nM, tBid 235 nM', ['E05']),
        ('Bax 0 nM, tBid 2350 nM', ['E06']),
        ('Bax 63 nM, tBid 0 nM', ['F01']),
        ('Bax 63 nM, tBid 0.24 nM', ['F02']),
        ('Bax 63 nM, tBid 2.35 nM', ['F03']),
        ('Bax 63 nM, tBid 23.5 nM', ['F04']),
        ('Bax 63 nM, tBid 235 nM', ['F05']),
        ('Bax 63 nM, tBid 2350 nM', ['F06']),
        ('Bax 252 nM, tBid 0 nM', ['G01']),
        ('Bax 252 nM, tBid 0.24 nM', ['G02']),
        ('Bax 252 nM, tBid 2.35 nM', ['G03']),
        ('Bax 252 nM, tBid 23.5 nM', ['G04']),
        ('Bax 252 nM, tBid 235 nM', ['G05']),
        ('Bax 252 nM, tBid 2350 nM', ['G06']),
        ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '130731_tBid_Bax_Timecourse.csv'))
triton_file = os.path.abspath(os.path.join(data_path,
                                        '130731_tBid_Bax_Triton.csv'))

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_wallac(timecourse_file)
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Post-Triton values
final_wells = read_wallac(triton_file)
final_wells = extract(wells_to_read, final_wells)
final_well_avgs = get_repeat_averages_by_well(final_wells)

# Initial timepoints
initial_wells_tc = get_first_points_by_well(timecourse_wells)

# Averages of raw timecourses across replicates
#(timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)

# Normalized timecourses
norm_wells = get_normalized_well_timecourses(
                    timecourse_wells, initial_wells_tc, final_well_avgs)

# Normalized timecourses, averaged
(norm_averages, norm_stds) = averages(norm_wells, layout)

# First timepoint shifted to 0 (better for fitting)
reset_norm_subset = reset_first_timepoint_to_zero(norm_averages)

# Pore timecourses
pores = get_average_pore_timecourses(norm_wells)

# pandas dataframe
#df = to_dataframe(norm_averages, norm_stds)

def plot_data():
    ion()

    figure()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    #figure()
    #plot_all(timecourse_averages, errors=timecourse_stds)
    #title("Raw timecourses, averaged")

    figure()
    plot_all(norm_wells, legend_position=0.7)
    title("Normalized timecourses")

    figure()
    plot_all(norm_averages, legend_position=0.7)
    title("Normalized timecourses, averaged")

    # Plot just the Bax 63 nM wells
    bax63_wells = ['Bax 63 nM, tBid 0 nM',
                   'Bax 63 nM, tBid 0.24 nM',
                   'Bax 63 nM, tBid 2.35 nM',
                   'Bax 63 nM, tBid 23.5 nM',
                   'Bax 63 nM, tBid 235 nM',
                   'Bax 63 nM, tBid 2350 nM',]
    figure()
    plot_all(extract(bax63_wells, norm_averages), legend_position=0.7)
    ylim((0, 1))
    title("63 nM Bax timecourses")

    # Plot just the Bax 252 nM wells
    bax252_wells = ['Bax 252 nM, tBid 0 nM',
                   'Bax 252 nM, tBid 0.24 nM',
                   'Bax 252 nM, tBid 2.35 nM',
                   'Bax 252 nM, tBid 23.5 nM',
                   'Bax 252 nM, tBid 235 nM',
                   'Bax 252 nM, tBid 2350 nM',]
    figure()
    plot_all(extract(bax252_wells, norm_averages), legend_position=0.7)
    ylim((0, 1))
    title("252 nM Bax timecourses")

    figure()
    plot_all(reset_norm_subset, legend_position=0.7)
    title("Norm., Avg., Reset to t = 0")

    figure()
    plot_all(pores)
    title("Avg. pores per liposome")

if __name__ == '__main__':
    plot_data()
