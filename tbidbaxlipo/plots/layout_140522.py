from tbidbaxlipo.util.plate_assay import *
import collections
import sys
import itertools
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
import numpy as np

layout = collections.OrderedDict([
        # ignore A12, a slight outlier
        ('Black seal, buffer',  ['A7', 'A8', 'A9', 'A10', 'A11']),
        ('Black seal, water',  ['B8', 'B9', 'B10', 'B11', 'B12']),
        ('No seal, buffer',  ['A1', 'A2', 'A3', 'A4', 'A5', 'A6']),
        ('No seal, water',  ['B1', 'B2', 'B3', 'B4', 'B5', 'B6']),
        ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '140522_ANTS_plt3881_background.txt'))

# Timecourse wells
timecourse_wells = read_flexstation_kinetics(timecourse_file)
"""The raw (unnormalized) timecourses."""

# Filter out only the wells that are in the layout dict
wells_to_read = list(itertools.chain(*layout.values()))
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Averages of raw timecourses across replicates
(timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)
"""Averages of raw timecourses."""

def plot_data():
    """Plots the data and various transformations of it."""

    # Timecourse wells
    figure()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    # Averages of raw timecourses across replicates
    figure()
    plot_all(timecourse_averages, errors=timecourse_stds)
    title("Raw timecourses, averaged")

    figure()
    for key, value in timecourse_averages.iteritems():
        plot(value[TIME], value[VALUE] / value[VALUE][0], label=key)
    title('Percentage change')
    legend(loc='center right')

if __name__ == '__main__':
    plt.ion()
    plot_data()
    sys.exit()

