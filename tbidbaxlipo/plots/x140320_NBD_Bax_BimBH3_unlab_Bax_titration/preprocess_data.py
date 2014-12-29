import sys
from os.path import dirname, abspath, join
import itertools
import collections
import tbidbaxlipo.data
from tbidbaxlipo.util.plate_assay import read_flexstation_kinetics, averages, \
                                         extract, subtract_background_set

layout = collections.OrderedDict([
        ('Bax 590 nM, NBD-Bax 96 nM',  ['A1', 'B1']),
        ('Bax 393 nM, NBD-Bax 96 nM',  ['A2', 'B2']),
        ('Bax 262 nM, NBD-Bax 96 nM',  ['A3', 'B3']),
        ('Bax 175 nM, NBD-Bax 96 nM',  ['A4', 'B4']),
        ('Bax 116 nM, NBD-Bax 96 nM',  ['A5', 'B5']),
        ('Bax 78 nM, NBD-Bax 96 nM',  ['A6', 'B6']),
        ('Bax 52 nM, NBD-Bax 96 nM',  ['A7', 'B7']),
        ('Bax 35 nM, NBD-Bax 96 nM',  ['A8', 'B8']),
        ('Bax 23 nM, NBD-Bax 96 nM',  ['A9', 'B9']),
        ('Bax 15 nM, NBD-Bax 96 nM',  ['A10', 'B10']),
        ('Bax 10 nM, NBD-Bax 96 nM',  ['A11', 'B11']),
        ('Bax 0 nM, NBD-Bax 96 nM',  ['A12', 'B12']),
        ('Bax 590 nM',  ['C1']),
        ('Bax 393 nM',  ['C2']),
        ('Bax 262 nM',  ['C3']),
        ('Bax 175 nM',  ['C4']),
        ('Bax 116 nM',  ['C5']),
        ('Bax 78 nM',  ['C6']),
        ('Bax 52 nM',  ['C7']),
        ('Bax 35 nM',  ['C8']),
        ('Bax 23 nM',  ['C9']),
        ('Bax 15 nM',  ['C10']),
        ('Bax 10 nM',  ['C11']),
        ('Bax 0 nM',  ['C12']),
    ])

data_path = dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = abspath(join(data_path,
                               '140320_NBD_Bax_BimBH3_unlab_Bax_titration.txt'))

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_flexstation_kinetics(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Averaged
(timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)

bax_bg_conditions = [
        'Bax 590 nM',
        'Bax 393 nM',
        'Bax 262 nM',
        'Bax 175 nM',
        'Bax 116 nM',
        'Bax 78 nM',
        'Bax 52 nM',
        'Bax 35 nM',
        'Bax 23 nM',
        'Bax 15 nM',
        'Bax 10 nM',
        'Bax 0 nM']

bax_bg_layout = extract(bax_bg_conditions, layout)
bax_bg_wells = extract([layout[cond][0] for cond in bax_bg_conditions],
                        timecourse_wells)

nbd_conditions = [
        'Bax 590 nM, NBD-Bax 96 nM',
        'Bax 393 nM, NBD-Bax 96 nM',
        'Bax 262 nM, NBD-Bax 96 nM',
        'Bax 175 nM, NBD-Bax 96 nM',
        'Bax 116 nM, NBD-Bax 96 nM',
        'Bax 78 nM, NBD-Bax 96 nM',
        'Bax 52 nM, NBD-Bax 96 nM',
        'Bax 35 nM, NBD-Bax 96 nM',
        'Bax 23 nM, NBD-Bax 96 nM',
        'Bax 15 nM, NBD-Bax 96 nM',
        'Bax 10 nM, NBD-Bax 96 nM',
        'Bax 0 nM, NBD-Bax 96 nM',]

nbd_layout = extract(nbd_conditions, layout)
nbd_wells = extract(nbd_conditions, timecourse_averages)
#nbd_wells = extract([layout[cond][0] for cond in nbd_conditions],
#                         timecourse_wells)

# Background subtracted
bgsub_wells = subtract_background_set(nbd_wells, bax_bg_wells)

# Background subtracted, averaged
#(bgsub_averages, bgsub_stds) = averages(bgsub_wells, layout)

# First timepoint shifted to 0 (better for fitting)
#reset_bgsub_means = reset_first_timepoint_to_zero(bgsub_norm_averages)
#Timecourses normalized, BG-subtracted, averaged, then with first point
#shifted to t = 0.
#reset_bgsub_sds = reset_first_timepoint_to_zero(bgsub_norm_stds)


