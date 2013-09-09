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
        #('Buffer', ['A01']),
        #('DMSO', ['A02']),
        #('365S', ['A03']),
        #('309S', ['A04']),
        #('310S', ['A05']),
        #('mDivi-1', ['A06']),
        ('tBid 47 nM', ['B01']),
        ('tBid 47 nM + DMSO', ['B02']),
        ('tBid 47 nM + analog B', ['B03']),
        ('tBid 47 nM + analog H', ['B04']),
        ('tBid 47 nM + 310S', ['B05']),
        ('tBid 47 nM + mDivi-1', ['B06']),
        ('Bax 126 nM', ['C01']),
        ('Bax 126 nM + DMSO', ['C02']),
        ('Bax 126 nM + analog B', ['C03']),
        ('Bax 126 nM + analog H', ['C04']),
        ('Bax 126 nM + 310S', ['C05']),
        ('Bax 126 nM + mDivi-1', ['C06']),
        ('tBid 47 nM + Bax 126 nM', ['D01', 'E01']),
        ('tBid 47 nM + Bax 126 nM + DMSO', ['D02', 'E02']),
        ('tBid 47 nM + Bax 126 nM + analog B', ['D03', 'E03']),
        ('tBid 47 nM + Bax 126 nM + analog H', ['D04', 'E04']),
        ('tBid 47 nM + Bax 126 nM + 310S', ['D05', 'E05']),
        ('tBid 47 nM + Bax 126 nM + mDivi-1', ['D06', 'E06']),
        ('Liposomes', ['F01']),
        ('DMSO', ['F02']),
        ('analog B', ['F03']),
        ('analog H', ['F04']),
        ('310S', ['F05']),
        ('mDivi-1', ['F06']),
        ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '130830_mDivi1_analogs_Timecourse.csv'))
triton_file = os.path.abspath(os.path.join(data_path,
                                        '130830_mDivi1_analogs_Triton.csv'))

ion()

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_wallac(timecourse_file)
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Initial timepoints
initial_wells_tc = get_first_points_by_well(timecourse_wells)

# Post-Triton values
final_wells = read_wallac(triton_file)
final_wells = extract(wells_to_read, final_wells)
final_well_avgs = get_repeat_averages_by_well(final_wells)

# Averages of raw timecourses across replicates
(timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)

# Normalized timecourses
norm_wells = get_normalized_well_timecourses(
        timecourse_wells, initial_wells_tc, final_well_avgs)

# Normalized timecourses, averaged
(norm_averages, norm_stds) = averages(norm_wells, layout)

# Pore timecourses
pores = get_average_pore_timecourses(norm_wells)

# First timepoint shifted to 0 (better for fitting)
reset_norm_subset = reset_first_timepoint_to_zero(norm_averages)

treatments_bax = extract(
        ['tBid 47 nM + Bax 126 nM',
        'tBid 47 nM + Bax 126 nM + DMSO',
        'tBid 47 nM + Bax 126 nM + analog B',
        'tBid 47 nM + Bax 126 nM + analog H',
        'tBid 47 nM + Bax 126 nM + 310S',
        'tBid 47 nM + Bax 126 nM + mDivi-1'],
        reset_norm_subset)

treatments_lipos_only = extract(
        ['Liposomes',
         'DMSO',
         'analog B',
         'analog H',
         '310S',
         'mDivi-1'],
        reset_norm_subset)

def plot_data():
    figure()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    figure(figsize=(10, 6))
    plot_all(timecourse_averages, errors=timecourse_stds, legend_position=0.6)
    title("Raw timecourses, averaged")

    figure()
    plot_all(norm_wells)
    title("Normalized timecourses")

    figure(figsize=(10, 6))
    plot_all(norm_averages, errors=norm_stds, legend_position=0.6)
    title("Normalized timecourses, averaged")

    figure()
    plot_all(pores)
    title("Avg. pores per liposome")

    figure(figsize=(10, 6))
    plot_all(reset_norm_subset, legend_position=0.6)
    title("Norm., Avg., Reset to t = 0")

    figure(figsize=(10, 6))
    plot_all(treatments_bax, legend_position=0.6)
    title('Dye release with tBid and Bax')

    figure(figsize=(10, 6))
    plot_all(treatments_lipos_only, legend_position=0.6)
    title('Dye release with liposomes only')

