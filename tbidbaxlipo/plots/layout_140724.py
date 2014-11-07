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
total_vol = 640. # uL
bax_stock = 190 # uL
bax_start_conc = ((bax_stock_conc * bax_stock)/total_vol) * 0.5

# Bid concs: 40, 10, 2.5, 0.625
bid06_labels = dose_series_replicate_list('Bid 0.63 nM, Bax', bax_start_conc, 2/3.,
                                         num_doses=12,
                                         start_row='C', end_row='C',
                                         start_col=1, end_col=12)
bid2_labels = dose_series_replicate_list('Bid 2.5 nM, Bax', bax_start_conc, 2/3.,
                                         num_doses=12,
                                         start_row='D', end_row='D',
                                         start_col=1, end_col=12)
bid10_labels = dose_series_replicate_list('Bid 10 nM, Bax', bax_start_conc, 2/3.,
                                         num_doses=12,
                                         start_row='E', end_row='E',
                                         start_col=1, end_col=12)
bid40_labels = dose_series_replicate_list('Bid 40 nM, Bax', bax_start_conc, 2/3.,
                                         num_doses=12,
                                         start_row='F', end_row='F',
                                         start_col=1, end_col=12)

bid06_layout = collections.OrderedDict(bid06_labels)
bid2_layout = collections.OrderedDict(bid2_labels)
bid10_layout = collections.OrderedDict(bid10_labels)
bid40_layout = collections.OrderedDict(bid40_labels)

timecourse_file = abspath(join(data_path, '140724_Bax_Bid_DPX_timecourse.txt'))
timecourse_wells = read_flexstation_kinetics(timecourse_file)
bid06_wells = extract(wells_from_layout(bid06_layout), timecourse_wells)
bid2_wells = extract(wells_from_layout(bid2_layout), timecourse_wells)
bid10_wells = extract(wells_from_layout(bid10_layout), timecourse_wells)
bid40_wells = extract(wells_from_layout(bid40_layout), timecourse_wells)

# Averages of raw timecourses across replicates
#(bax_averages, bax_stds) = averages(timecourse_wells, bax_layout)

# In this experiment, we first prepare liposomes at the same final
# concentration as our wells that will contain Bid and Bax; lyse them with
# Triton; then sequentially add back DPX and measure the quenching at each
# DPX concentration.  Each addition of DPX is read and hence goes in a
# separate file. This list contains the names of files containing the data
# for the standard curve:
dpx_std_file_list = [
        '140724_Bax_Bid_DPX_1_0uL.txt',
        '140724_Bax_Bid_DPX_2_1uL.txt',
        '140724_Bax_Bid_DPX_3_2uL.txt',
        '140724_Bax_Bid_DPX_4_4uL.txt',
        '140724_Bax_Bid_DPX_5_6uL.txt',
        '140724_Bax_Bid_DPX_6_8uL.txt',
        '140724_Bax_Bid_DPX_7_10uL.txt']
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
bid06_requench_wells = [row + col for (col, row) in itertools.product(
                            [str(i) for i in range(1, 13)],
                            ['C'])]
bid2_requench_wells = [row + col for (col, row) in itertools.product(
                            [str(i) for i in range(1, 13)],
                            ['D'])]
bid10_requench_wells = [row + col for (col, row) in itertools.product(
                            [str(i) for i in range(1, 13)],
                            ['E'])]
bid40_requench_wells = [row + col for (col, row) in itertools.product(
                            [str(i) for i in range(1, 13)],
                            ['F'])]

fmax_filename = abspath(join(data_path, '140724_Bax_Bid_DPX_triton.txt'))

def main():
    plt.close('all')

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

    # If we plot the BG values as a function of DPX concentration, we can
    # identify the outliers (probably a result of a piece of lint/debris on the
    # surface of the plate)
    plt.figure()
    plt.plot(dpx_concs, bg_matrix)
    plt.title('DPX background fluorescence')
    plt.xlabel('DPX concentration')
    plt.ylabel('ANTS (RFU)')

    # We set the outlier values to NaN
    bg_matrix[bg_matrix > 70] = np.nan

    # Now, we take the mean at each DPX concentration. These are the values we
    # will subtract from every well before calculating quenching.
    bg_avgs = np.nanmean(bg_matrix, axis=1)

    # We replot, with means and error bars
    plt.figure()
    plt.errorbar(dpx_concs, np.nanmean(bg_matrix, axis=1),
                 yerr=np.nanstd(bg_matrix, axis=1), color='k', linewidth=2)
    plt.title('DPX background fluorescence, average')
    plt.xlabel('DPX concentration')
    plt.ylabel('ANTS (RFU)')

    (i_avgs, i_sds, fmax_avg, fmax_sd) = \
            quenching_std_curve(dpx_std_file_list, dpx_std_wells, dpx_concs,
                                bg_avgs=bg_avgs)
    (ka, kd) = fit_std_curve(i_avgs, i_sds, dpx_concs)
    qd = get_quenching_dict(i_avgs, i_sds, dpx_vols_added)
    final_q = qd[dpx_vols_added[-1]]
    q_outs = np.array(qd.values())

    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bid06_requench_wells,
                                         final_q, final_bg=bg_avgs[-1])
    requenching_analysis(dpx_std_file_list, bid06_requench_wells, dpx_concs,
            q_outs, fmax_avgs, fmax_sds, None, None, None, bg_avgs)

    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bid40_requench_wells,
                                         final_q, final_bg=bg_avgs[-1])
    requenching_analysis(dpx_std_file_list, bid40_requench_wells, dpx_concs,
            q_outs, fmax_avgs, fmax_sds, None, None, None, bg_avgs)



if __name__ == '__main__':
    main()
