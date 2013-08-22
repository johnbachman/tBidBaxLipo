from tbidbaxlipo.util.plate_assay import *
import collections
import sys
import os
import tbidbaxlipo.data

# In the course of doing this experiment, I learned that trying to pipette 50
# uL into the half-well area plate leads to spillage and error. 20 uL had
# no sign of spilling.

well_names = ['%s%.2d' % (row, col)
              for row in ['C','D','E','F','G','H']
              for col in range(1,13)]

layout = collections.OrderedDict([('Dye 20 uL', well_names)])

ion()

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '130822_wallac_dispenser_20ul.csv'))

timecourse_wells = read_wallac(timecourse_file)
"""The raw (unnormalized) timecourses."""

reset_wells = reset_well_times(timecourse_wells)
"""Reset each kinetic timecourse to have t = 0."""

trunc_wells = truncate_timecourses(reset_wells, 15)
"""Truncate the first 15 (of 30) points."""

# Calculate the mean, SD, and CV for the timecourse in each well
means = []
sds = []
cvs = []
for well in trunc_wells.keys():
    means.append(np.mean(trunc_wells[well][VALUE]))
    sds.append(np.std(trunc_wells[well][VALUE]))
for i, mean in enumerate(means):
    cvs.append(100 * (sds[i] / float(mean)))

# Now, find the pipetting error
global_mean = np.mean(means)
global_sd = np.std(means)
global_cv = 100 * (global_sd / global_mean)

def main():
    figure()
    plot_all(reset_wells)
    title("20uL timecourses, raw")

    figure()
    plot_all(trunc_wells)
    title("20ul timecourses, first 15 points truncated")

    figure()
    hist(means, color='r')
    title("Distribution of means for individual wells")
    xlabel("Timecourse mean")
    ylabel("Count")

    figure()
    hist(sds, color='r')
    title("Distribution of SDs for individual wells")
    xlabel("Timecourse standard deviation")
    ylabel("Count")

    figure()
    hist(cvs, color='r')
    title("Distribution of CVs for individual wells")
    xlabel("Timecourse coeff. of variation")
    ylabel("Count")

def print_statistics():
    print "Mean across dispensed wells: %f" % global_mean
    print "SD across dispensed wells:   %f" % global_sd
    print "CV across dispensed wells:   %f" % global_cv

def plot_edge_effects():
    # Now we plot across rows to look for edge effects
    wells_by_row = []
    for row in ['C','D','E','F','G','H']:
        row_wells = []
        for col in range(1,13):
            well_name = '%s%.2d' % (row, col)
            row_wells.append(well_name)
        wells_by_row.append(row_wells)

    figure()
    for row in wells_by_row:
        row_vals = [np.mean(trunc_wells[well_name][VALUE])
                    for well_name in row]
        plot(row_vals)


