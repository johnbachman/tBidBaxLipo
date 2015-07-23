from tbidbaxlipo.util.plate_assay import \
        extract, read_flexstation_kinetics, subtract_background_set, averages, \
        TIME, VALUE
import collections
import tbidbaxlipo.data
import sys
from os.path import dirname, abspath, join
import itertools
import numpy as np

layout = collections.OrderedDict([
        ('Bax 185 nM, Lipos 1 mg/ml',  ['A1']),
        ('Bax 185 nM, Lipos 0.5 mg/ml',  ['A2']),
        ('Bax 185 nM, Lipos 0.25 mg/ml',  ['A3']),
        ('Bax 185 nM, Lipos 0.125 mg/ml',  ['A4']),
        ('Bax 185 nM, Lipos 0.063 mg/ml',  ['A5']),
        ('Bax 185 nM, Lipos 0.031 mg/ml',  ['A6']),
        ('Bax 185 nM, Lipos 0.016 mg/ml',  ['A7']),
        ('Bax 185 nM, Lipos 0.008 mg/ml',  ['A8']),
        ('Bax 185 nM, Lipos 0.004 mg/ml',  ['A9']),
        ('Bax 185 nM, Lipos 0.002 mg/ml',  ['A10']),
        ('Bax 185 nM, Lipos 0.001 mg/ml',  ['A11']),
        ('Bax 185 nM, Lipos 0 mg/ml',  ['A12']),
        ('Lipos 1 mg/ml',  ['B1']),
        ('Lipos 0.5 mg/ml',  ['B2']),
        ('Lipos 0.25 mg/ml',  ['B3']),
        ('Lipos 0.125 mg/ml',  ['B4']),
        ('Lipos 0.063 mg/ml',  ['B5']),
        ('Lipos 0.031 mg/ml',  ['B6']),
        ('Lipos 0.016 mg/ml',  ['B7']),
        ('Lipos 0.008 mg/ml',  ['B8']),
        ('Lipos 0.004 mg/ml',  ['B9']),
        ('Lipos 0.002 mg/ml',  ['B10']),
        ('Lipos 0.001 mg/ml',  ['B11']),
        ('Lipos 0 mg/ml',  ['B12']),
    ])
data_path = dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = abspath(join(data_path,
                              '140318_NBD_Bax_BimBH3_lipo_titration.txt'))

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_flexstation_kinetics(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = extract(wells_to_read, timecourse_wells)

lipo_bg_conditions = [
        'Lipos 1 mg/ml',
        'Lipos 0.5 mg/ml',
        'Lipos 0.25 mg/ml',
        'Lipos 0.125 mg/ml',
        'Lipos 0.063 mg/ml',
        'Lipos 0.031 mg/ml',
        'Lipos 0.016 mg/ml',
        'Lipos 0.008 mg/ml',
        'Lipos 0.004 mg/ml',
        'Lipos 0.002 mg/ml',
        'Lipos 0.001 mg/ml',
        'Lipos 0 mg/ml']
lipo_bg_layout = extract(lipo_bg_conditions, layout)
lipo_bg_wells = extract([layout[cond][0] for cond in lipo_bg_conditions],
                        timecourse_wells)

bax_lipo_conditions = [
        'Bax 185 nM, Lipos 1 mg/ml',
        'Bax 185 nM, Lipos 0.5 mg/ml',
        'Bax 185 nM, Lipos 0.25 mg/ml',
        'Bax 185 nM, Lipos 0.125 mg/ml',
        'Bax 185 nM, Lipos 0.063 mg/ml',
        'Bax 185 nM, Lipos 0.031 mg/ml',
        'Bax 185 nM, Lipos 0.016 mg/ml',
        'Bax 185 nM, Lipos 0.008 mg/ml',
        'Bax 185 nM, Lipos 0.004 mg/ml',
        'Bax 185 nM, Lipos 0.002 mg/ml',
        'Bax 185 nM, Lipos 0.001 mg/ml',
        'Bax 185 nM, Lipos 0 mg/ml', ]
bax_lipo_layout = extract(bax_lipo_conditions, layout)
bax_lipo_wells = extract([layout[cond][0] for cond in bax_lipo_conditions],
                         timecourse_wells)

# Normalized and background subtracted
bgsub_wells = subtract_background_set(bax_lipo_wells, lipo_bg_wells)

(bgsub_averages, bgsub_sds) = averages(bgsub_wells, bax_lipo_layout)

# First timepoint shifted to 0 (better for fitting)
#reset_bgsub_means = reset_first_timepoint_to_zero(bgsub_norm_averages)
"""Timecourses normalized, BG-subtracted, averaged, then with first point
shifted to t = 0."""
#reset_bgsub_sds = reset_first_timepoint_to_zero(bgsub_norm_stds)

lipo_conc_conv_factor = 15.502 # 1 mg/ml ~= 15.502 nM liposomes
bg_tc = bgsub_averages['Bax 185 nM, Lipos 0 mg/ml'][VALUE]
bg_time = bgsub_averages['Bax 185 nM, Lipos 0 mg/ml'][TIME]
data_to_fit = []
lipo_concs_to_fit = []
lipo_mgs_to_fit = [1., 0.5, 0.25, 0.125, 0.063, 0.031, 0.016, 0.008, 0]
for conc_name in bgsub_averages.keys():
    lipo_mg = float(conc_name.split()[4])
    if not lipo_mg in lipo_mgs_to_fit:
        continue
    lipo_concs_to_fit.append(lipo_mg * lipo_conc_conv_factor)
    t = bgsub_averages[conc_name][TIME]
    v = bgsub_averages[conc_name][VALUE]
    v_bg = v / bg_tc
    data_to_fit.append(v_bg)
lipo_concs_to_fit = np.array(lipo_concs_to_fit)

# Variable containing estimate of experimental error
data_sigma = [0.1]
