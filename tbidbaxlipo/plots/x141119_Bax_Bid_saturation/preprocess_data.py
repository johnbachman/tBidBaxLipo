import sys
from os.path import dirname, abspath, join
import collections
import numpy as np
import tbidbaxlipo.data
from tbidbaxlipo.util.plate_assay import read_flexstation_kinetics, \
                        add_offset_vector, row_wells, extract, averages, \
                        subtract_background, TIME, VALUE
from scipy.stats import linregress

# All wells with liposomes ~0.15 mg/mL

# The path to the CSV data file from the plate reader
data_path = dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = abspath(join(data_path,
                               '141119_Bax_Bid_saturation_timecourse.txt'))

# The concentrations of Bax, including NBD-C126-Bax and competitor.
# These are given in the order of the wells, e.g., column 1 has 1000 nM
# and column 11 has 0 nM.
bax_concs = np.array([1000., 500., 250., 125., 62.5, 31.2, 15.6, 7.8,
                      3.9, 2.0, 0.]) + 25.
# Concentrations of cBid.
bid_concs = np.array([0.0, 2.5, 5., 10., 20., 40., 80.])

# Parse the CSV file from the instrument into a dict of timecourses,
# indexed by well name.
timecourse_wells = read_flexstation_kinetics(timecourse_file)

def time_offset_vector():
    """Calculate the time offset vector for this experiment.

    Because Bax was manually added to the wells and mixed row by row, the
    actual start of the experiment was slightly different for each row of
    wells.  During the experiment the mixing time for each row was noted and is
    here used to calculate the offset for each row of wells. This way the wells
    that have been incubating longer are recorded as starting at an earlier
    time than the wells that were mixed immediately before reading.
    """

    # The time, in seconds, from the mixing of the first row to the first
    # read of the plate
    read_time = 320 # seconds (5:20)
    # Elapsed time of mixing since first mixed row
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
    # In other words, the first timepoint (read time t = 0) for row A is at
    # (320 - 0) = 320 seconds after mixing
    # whereas the first timepoint for row H is at
    # (320 - 150) = 170 seconds after mixing
    for row in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
        for col in range(1, 13):
            well = '%s%s' % (row, col)
            offset_dict[well] = read_time - row_times[row]
    return offset_dict

# Add the time offset to the time coordinates for the data
timecourse_wells = add_offset_vector(timecourse_wells, time_offset_vector())

# No NBD or lipos (buffer only)
# The buffer only wells are used to establish the background fluorescence
# of the buffer/well.
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
# Produces average timecourses across replicates for the three background
# conditions
(bg_averages, bg_std) = averages(timecourse_wells, bg_layout)

# To do the correct background subtraction we need to account for the
# fluorescence of both the buffer itself and the liposomes. However, I
# unfortunately neglected to include a control condition with liposomes and no
# NBD-Bax. Fortunately, I can back-calculate the fluorescence of the liposomes
# alone by sutracting the NBD-Bax + buffer condition from the NBD-Bax + lipos +
# buffer condition. This value can then be added back to the buffer only
# condition to get the background vector (buffer + lipos) that will be
# subtracted from all NBD-Bax timecourses.
#
# However, note that since the NBD-Bax + lipos + buffer condition will include
# fluorescence increases due to spontaneous activation and insertion of Bax,
# whereas the NBD-Bax + buffer (no liposomes) will not incorporate this
# increase. Therefore ONLY the initial timepoints (say, the average of the
# first three timepoints) should be used to establish the liposome-only
# background.

# The difference between NBD-Bax and background (buffer only) gives fluorescence
# of NBD-Bax alone, with no liposomes (we don't use this value, but it can
# be inspected for reference).
bg_diff = {}
bg_diff['NBD-Bax only'] = []
bg_diff['NBD-Bax only'].append(bg_averages['No lipos'][TIME])
bg_diff['NBD-Bax only'].append(bg_averages['No lipos'][VALUE] -
                          bg_averages['No NBD or lipos'][VALUE])
# Calculate the average fluorescence timecourse of
# (NBD-Bax + lipos + buffer) - (NBD-Bax + buffer) = lipos
bg_diff['Liposomes only'] = []
bg_diff['Liposomes only'].append(bg_averages['NBD and lipos'][TIME])
bg_diff['Liposomes only'].append(bg_averages['NBD and lipos'][VALUE] -
                           bg_averages['No lipos'][VALUE])

# Calculate the average initial liposome-only fluorescence (before NBD-Bax
# starts to increase due to spontaneous insertion).
# Liposomes alone has fluorescence of roughly 2.54.
lipo_bg_offset = np.mean(bg_diff['Liposomes only'][VALUE][0:3])

# Add the calculated lipo only background to the buffer-only background vector
bg_diff['Est. buffer + lipos'] = []
bg_diff['Est. buffer + lipos'].append(bg_averages['No NBD or lipos'][TIME])
bg_diff['Est. buffer + lipos'].append(bg_averages['No NBD or lipos'][VALUE] + \
                                      lipo_bg_offset)

# So need to subtract background and 2.54 from each curve to get
# baseline fluorescence. The other curves are normalized against this.
bg_to_subtract = bg_averages['No NBD or lipos'][VALUE] + lipo_bg_offset

# Subtract background from all wells
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

#data_matrix = np.zeros((len(bid_concs), len(bax_concs), len(bg_averages))

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
time = t
