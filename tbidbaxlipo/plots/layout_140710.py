from tbidbaxlipo.util.dpx_assay import *
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import tbidbaxlipo.data
from tbidbaxlipo.util.plate_assay import *
from os.path import abspath, join
import collections

plt.ion()

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)

# Calculate starting dose for first Bax dilution series
bax_stock_conc = 5900. # nM
total_vol = 170. # uL
bax_stock1 = 46.1 # uL
bax_stock2 = 38.4 # uL
bax_start_conc1 = ((bax_stock_conc * bax_stock1)/total_vol) * 0.5
bax_start_conc2 = ((bax_stock_conc * bax_stock2)/total_vol) * 0.5


bax_labels1 = dose_series_replicate_list('Bax', bax_start_conc1, 2/3.,
                                         num_doses=12,
                                         start_row='C', end_row='C',
                                         start_col=1, end_col=12)
bax_labels2 = dose_series_replicate_list('Bax', bax_start_conc2, 2/3.,
                                         num_doses=12,
                                         start_row='D', end_row='D',
                                         start_col=1, end_col=12)

bax_layout = collections.OrderedDict(bax_labels1 + bax_labels2)

cec_layout = collections.OrderedDict(dose_series_replicate_list(
                    'CecropinA', 20000, 2/3., num_doses=12,
                    start_row='E', end_row='E', start_col=1, end_col=12))

timecourse_file = abspath(join(data_path, '140710_Bax_43C_DPX_timecourse.txt'))
timecourse_wells = read_flexstation_kinetics(timecourse_file)
bax_wells = extract(wells_from_layout(bax_layout), timecourse_wells)
cec_wells = extract(wells_from_layout(cec_layout), timecourse_wells)

# Averages of raw timecourses across replicates
(bax_averages, bax_stds) = averages(timecourse_wells, bax_layout)
(cec_averages, cec_stds) = averages(timecourse_wells, cec_layout)

# In this experiment, we first prepare liposomes at the same final
# concentration as our wells that will contain Bid and Bax; lyse them with
# Triton; then sequentially add back DPX and measure the quenching at each
# DPX concentration.  Each addition of DPX is read and hence goes in a
# separate file. This list contains the names of files containing the data
# for the standard curve:
dpx_std_file_list = [
        '140710_Bax_43C_DPX_1_0uL_fixed.txt',
        '140710_Bax_43C_DPX_2_1uL.txt',
        '140710_Bax_43C_DPX_3_2uL.txt',
        '140710_Bax_43C_DPX_4_4uL.txt',
        '140710_Bax_43C_DPX_5_6uL.txt',
        '140710_Bax_43C_DPX_6_8uL.txt',
        '140710_Bax_43C_DPX_7_10uL.txt']
dpx_std_file_list = [abspath(join(data_path, filename))
                     for filename in dpx_std_file_list]


# These are the wells in the plate that contain the liposome-only solution
# used to calculate the standard curve:
dpx_std_wells = ['B%d' % well_num for well_num in range(1, 13)]

# The list of DPX volumes added at each step
dpx_vol_steps = [0., 1., 2., 4., 6., 8., 10.]

# The cumulative volumes of stock DPX added at each step
dpx_vols_added = np.cumsum(dpx_vol_steps)

# We know how much of the stock DPX we added at each quenching step, but we
# have to calculate the actual concentrations:
dpx_concs = calc_dpx_concs(dpx_vol_steps, starting_well_vol=100.)

# The wells containing the Bid/Bax treated liposomes
bax_requench_wells = [row + col for (col, row) in itertools.product(
                            [str(i) for i in range(1, 13)],
                            ['C', 'D'])]
# The wells containing the Cecropin A treated liposomes
cec_requench_wells = [row + col for (col, row) in itertools.product(
                            [str(i) for i in range(1, 13)],
                            ['E'])]


fmax_filename = abspath(join(data_path, '140710_Bax_43C_DPX_triton.txt'))

## DPX BACKGROUND FLUORESCENCE

# These are the wells in the plate that contain buffer only but to which
# DPX is also added; used for background subtraction
bg_wells = ['A%d' % well_num for well_num in range(1, 13)]

# Iterate over the DPX files, collecting bg values
bg_matrix = np.zeros((len(dpx_std_file_list), len(bg_wells)))
for file_index, file in enumerate(dpx_std_file_list):
    wells = read_flexstation_kinetics(file)
    for well_index, bg_well in enumerate(bg_wells):
        bg_matrix[file_index, well_index] = np.mean(wells[bg_well][VALUE])

# If we plot the BG values as a function of DPX concentration, we can identify
# the outliers (probably a result of a piece of lint/debris on the surface of
# the plate)
plt.figure()
plt.plot(dpx_concs, bg_matrix)
plt.title('DPX background fluorescence')
plt.xlabel('DPX concentration')
plt.ylabel('ANTS (RFU)')

# We set the outlier values to NaN
bg_matrix[bg_matrix > 70] = np.nan

# Now, we take the mean at each DPX concentration. These are the values we will
# subtract from every well before calculating quenching.
bg_avgs = np.nanmean(bg_matrix, axis=1)

# We replot, with means and error bars
plt.figure()
plt.errorbar(dpx_concs, np.nanmean(bg_matrix, axis=1),
             yerr=np.nanstd(bg_matrix, axis=1), color='k', linewidth=2)
plt.title('DPX background fluorescence')
plt.xlabel('DPX concentration')
plt.ylabel('ANTS (RFU)')

## DILUTION

# How much is the signal affected by the simple dilution of the wells?
# To address this, added a row (F) to the plate that had liposomes lysed with
# Triton, to which buffer (no DPX) was added at the same volumes as the DPX.

dilution_wells = ['F%d' % well_num for well_num in range(1, 13)]

# Iterate over the DPX files, collecting bg values
dilution_matrix = np.zeros((len(dpx_std_file_list), len(dilution_wells)))
for file_index, file in enumerate(dpx_std_file_list):
    wells = read_flexstation_kinetics(file)
    for well_index, d_well in enumerate(dilution_wells):
        dilution_matrix[file_index, well_index] = np.mean(wells[d_well][VALUE])

# Normalize each value on a well-by-well basis according to the value with
# no dilution.
# The no dilution wells are in the first row (first file):
no_dilution_wells = dilution_matrix[0, :]
# Divide all wells by the values in the no dilution wells
dilution_matrix = dilution_matrix / no_dilution_wells

# Look for outliers
plt.figure()
plt.plot(dpx_vols_added, dilution_matrix)
plt.title('Effect of dilution on fluorescence')
plt.xlabel('Volume added')
plt.ylabel('ANTS (RFU)')

# No real outliers here.
# Plot means with error bars:
plt.figure()
plt.errorbar(dpx_vols_added, np.mean(dilution_matrix, axis=1),
             yerr=np.std(dilution_matrix, axis=1), color='k', linewidth=2)
plt.title('Effect of dilution on fluorescence')
plt.xlabel('Volume added')
plt.ylabel('ANTS (RFU)')

# Note that dilution curve first goes up slightly before it goes down.
# Since it is unlikely that this is due to the 1 uL of buffer that is added,
# it seems likely that since the liposomes were lysed immediately before the
# plate was read, the Triton hadn't yet finished lysing all of the liposomes.

# I'm also not sure it's necessary to break up the dilution effect and the
# quenching effect into separate components?

# For now I will proceed just by doing the background subtraction, and will set
# aside the dilution issue.

# Background subtraction: subtract the background for that DPX concentration
# from each well, so that the baseline fluorescence will be due only to the
# unquenched ANTS fluorescence.

# Need to look for outliers in the standard curve!




# But the question remains, how much of the 

#plot_endpoints_vs_dose(bax_averages, bax_layout)
#plot_endpoints_vs_dose(cec_averages, cec_layout)

(i_avgs, i_sds, fmax_avg, fmax_sd) = \
        quenching_std_curve(dpx_std_file_list,
                                    dpx_std_wells, dpx_concs, bg_avgs=bg_avgs)
(ka, kd) = fit_std_curve(i_avgs, i_sds, dpx_concs)

qd = get_quenching_dict(i_avgs, i_sds, dpx_vols_added)

final_q = qd[dpx_vols_added[-1]]

(fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bax_requench_wells,
                                     final_q, final_bg=bg_avgs[-1])

q_outs = np.array(qd.values())

requenching_analysis(dpx_std_file_list, bax_requench_wells, dpx_concs, q_outs,
                     fmax_avgs, fmax_sds, None, None, None, bg_avgs)
