import sys
from os.path import dirname, abspath, join
import collections
import numpy as np
import tbidbaxlipo.data
from tbidbaxlipo.util.plate_assay import read_flexstation_kinetics, \
                        add_offset_vector, row_wells, extract, averages, \
                        subtract_background, TIME, VALUE

# All wells with liposomes ~0.15 mg/mL

data_path = dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = abspath(join(data_path,
                               '141119_Bax_Bid_saturation_timecourse.txt'))

bax_concs = np.array([1000., 500., 250., 125., 62.5, 31.2, 15.6, 7.8,
                      3.9, 2.0, 0.]) + 25.


# Timecourse wells
timecourse_wells = read_flexstation_kinetics(timecourse_file)
"""The raw (unnormalized) timecourses."""

# Specify the offset vector for this experiment
def time_offset_vector():
    read_time = 320 # seconds (5:20)
    row_times = {'A': 0,
                 'B': 20,
                 'C': 37,
                 'D': 60,
                 'E': 80,
                 'F': 100,
                 'G': 130,
                 'H': 150,
             }
    offset_dict = {}
    # New time vector is TIME + offset, where
    # offset = (read_time - row_time)
    for row in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
        for col in range(1, 13):
            well = '%s%s' % (row, col)
            offset_dict[well] = read_time - row_times[row]
    return offset_dict

# Add the time offset to the time coordinates for the data
timecourse_wells = add_offset_vector(timecourse_wells, time_offset_vector())

# No NBD or lipos (buffer only)
no_nbd_lipo_well_names = ['B12', 'D12', 'F12'] # Ignoring H12 as an outlier
no_nbd_lipo_wells = extract(no_nbd_lipo_well_names, timecourse_wells)

# NBD-Bax but no lipos
no_lipo_well_names = ['A12', 'C12', 'E12', 'G12']
no_lipo_wells = extract(no_lipo_well_names, timecourse_wells)

bg_layout = collections.OrderedDict([
         ('No lipos', ['A12', 'C12', 'E12', 'G12']),
         ('No NBD or lipos', ['B12', 'D12', 'F12']),
         ('NBD and lipos', ['G9', 'G10', 'G11']),
        ])
# Produces average timecourses for the three BG conditions separately
(bg_averages, bg_std) = averages(timecourse_wells, bg_layout)

# Difference between NBD-Bax and background (buffer only) gives fluorescence
# of NBD-Bax alone
bg_diff = {}
bg_diff['BG diff'] = []
bg_diff['BG diff'].append(bg_averages['No lipos'][TIME])
bg_diff['BG diff'].append(bg_averages['No lipos'][VALUE] -
                          bg_averages['No NBD or lipos'][VALUE])
bg_diff['BG diff2'] = []
bg_diff['BG diff2'].append(bg_averages['NBD and lipos'][TIME])
bg_diff['BG diff2'].append(bg_averages['NBD and lipos'][VALUE] -
                           bg_averages['No lipos'][VALUE])
bg_diff['BG diff3'] = []
bg_diff['BG diff3'].append(bg_averages['No NBD or lipos'][TIME])
bg_diff['BG diff3'].append(bg_averages['No NBD or lipos'][VALUE] + 2.5)

# So liposomes alone has fluorescence of roughly 2.55
# So need to subtract background and 2.55 from each curve to get
# baseline fluorescence. The other curves are normalized against this.

bg_to_subtract = bg_averages['No NBD or lipos'][VALUE] + 2.5

bg_sub_wells = subtract_background(timecourse_wells, bg_to_subtract)

# These dicts allow us to plot and analyze the different Bid conditions
# separately:
bid_80 = extract(row_wells('A', 11), bg_sub_wells)
bid_40 = extract(row_wells('B', 11), bg_sub_wells)
bid_20 = extract(row_wells('C', 11), bg_sub_wells)
bid_10 = extract(row_wells('D', 11), bg_sub_wells)
bid_5 = extract(row_wells('E', 11), bg_sub_wells)
bid_2 = extract(row_wells('F', 11), bg_sub_wells)
bid_0 = extract(row_wells('G', 11), bg_sub_wells)
bim_bh3 = extract(row_wells('H', 11), bg_sub_wells)

# Fit the  20 nM Bid condition
data_to_fit = []
numpts = 20
for well in bid_20.keys():
    t = bid_20[well][TIME]
    v = bid_20[well][VALUE]
    lin_fit = linregress(t[:numpts], v[:numpts])
    intercept = lin_fit[1]
    norm_data = (v / intercept)
    data_to_fit.append(norm_data)
