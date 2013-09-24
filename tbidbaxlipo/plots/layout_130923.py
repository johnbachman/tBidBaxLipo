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
        ('Liposomes + DMSO', ['A01']),
        ('Liposomes + mdivi-1', ['A02']),
        ('Liposomes + 309S', ['A03']),
        ('Liposomes + 310S', ['A04']),
        ('Liposomes + 365S', ['A05']),
        ('cBid 47 nM + DMSO', ['B01']),
        ('cBid 47 nM + mdivi-1', ['B02']),
        ('cBid 47 nM + 309S', ['B03']),
        ('cBid 47 nM + 310S', ['B04']),
        ('cBid 47 nM + 365S', ['B05']),
        ('Bax 126 nM + DMSO', ['C01']),
        ('Bax 126 nM + mdivi-1', ['C02']),
        ('Bax 126 nM + 309S', ['C03']),
        ('Bax 126 nM + 310S', ['C04']),
        ('Bax 126 nM + 365S', ['C05']),
        ('cBid 47 nM + Bax 126 nM + DMSO', ['D01']),
        ('cBid 47 nM + Bax 126 nM + mdivi-1', ['D02']),
        ('cBid 47 nM + Bax 126 nM + 309S', ['D03']),
        ('cBid 47 nM + Bax 126 nM + 310S', ['D04']),
        ('cBid 47 nM + Bax 126 nM + 365S', ['D05']),
    ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                 '130923_fitc_dextran_mdivi1_timecourse.csv'))
triton_file = os.path.abspath(os.path.join(data_path,
                                 '130923_fitc_dextran_mdivi1_triton.csv'))

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

# Normalized timecourses
norm_wells = get_normalized_well_timecourses(
        timecourse_wells, initial_wells_tc, final_well_avgs)

# Normalized timecourses, averaged
(norm_averages, norm_stds) = averages(norm_wells, layout)

# Pore timecourses
pores = get_average_pore_timecourses(norm_wells)

# First timepoint shifted to 0 (better for fitting)
reset_norm_subset = reset_first_timepoint_to_zero(norm_averages)

treatments_cbid_bax = extract(
        ['cBid 47 nM + Bax 126 nM + DMSO',
        'cBid 47 nM + Bax 126 nM + mdivi-1',
        'cBid 47 nM + Bax 126 nM + 309S',
        'cBid 47 nM + Bax 126 nM + 310S',
        'cBid 47 nM + Bax 126 nM + 365S'],
        norm_averages)

dmso = extract(
        ['Liposomes + DMSO',
         'cBid 47 nM + DMSO',
         'Bax 126 nM + DMSO',
         'cBid 47 nM + Bax 126 nM + DMSO'],
        norm_averages)

mdivi1 = extract(
        ['Liposomes + mdivi-1',
         'cBid 47 nM + mdivi-1',
         'Bax 126 nM + mdivi-1',
         'cBid 47 nM + Bax 126 nM + mdivi-1'],
        norm_averages)

m309S = extract(
        ['Liposomes + 309S',
         'cBid 47 nM + 309S',
         'Bax 126 nM + 309S',
         'cBid 47 nM + Bax 126 nM + 309S'],
        norm_averages)

m310S = extract(
        ['Liposomes + 310S',
         'cBid 47 nM + 310S',
         'Bax 126 nM + 310S',
         'cBid 47 nM + Bax 126 nM + 310S'],
        norm_averages)

m365S = extract(
        ['Liposomes + 365S',
         'cBid 47 nM + 365S',
         'Bax 126 nM + 365S',
         'cBid 47 nM + Bax 126 nM + 365S'],
        norm_averages)

def plot_data():
    def myfig():
        figure(figsize=(10, 6))

    myfig()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    figure(figsize=(10,6))
    plot_all(norm_averages, legend_position=0.7)
    title("Normalized timecourses, averaged")

    myfig()
    plot_all(pores, legend_position=0.7)
    title("Avg. pores per liposome")

    myfig()
    plot_all(treatments_cbid_bax, legend_position=0.7)
    title('Treatment comparison, cBid + Bax')

    myfig()
    plot_all(dmso, legend_position=0.7)
    title('Dye release, DMSO')

    myfig()
    plot_all(mdivi1, legend_position=0.7)
    title('Dye release, mdivi-1')

    myfig()
    plot_all(m309S, legend_position=0.7)
    title('Dye release, 309S')

    myfig()
    plot_all(m310S, legend_position=0.7)
    title('Dye release, 310S')

    myfig()
    plot_all(m365S, legend_position=0.7)
    title('Dye release, 365S')

if __name__ == '__main__':
    plot_data()
