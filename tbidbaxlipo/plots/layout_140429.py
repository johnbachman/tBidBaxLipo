from tbidbaxlipo.util.dpx_assay import *
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import tbidbaxlipo.data
from tbidbaxlipo.util.plate_assay import *
from os.path import abspath, join
import collections

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)

bax_labels1 = dose_series_replicate_list('Bax', 900., 2/3., num_doses=12,
                                        start_row='D', end_row='D',
                                        start_col=1, end_col=12)
bax_labels2 = dose_series_replicate_list('Bax', 750., 2/3., num_doses=12,
                                        start_row='E', end_row='E',
                                        start_col=1, end_col=12)
bax_layout = collections.OrderedDict(bax_labels1 + bax_labels2)

cec_layout = collections.OrderedDict(dose_series_replicate_list(
                    'CecropinA', 50000, 2/3., num_doses=12,
                    start_row='C', end_row='C', start_col=1, end_col=12))

timecourse_file = abspath(join(data_path, '140429_Bax_DPX_timecourse.txt'))
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
        '140429_Bax_DPX_post_triton.txt',
        '140429_Bax_DPX_1_1uL.txt',
        '140429_Bax_DPX_2_2uL.txt',
        '140429_Bax_DPX_3_4uL.txt',
        '140429_Bax_DPX_4_6uL.txt',
        '140429_Bax_DPX_5_8uL.txt',
        '140429_Bax_DPX_6_10uL.txt',
        '140429_Bax_DPX_7_10uL.txt']
dpx_std_file_list = [abspath(join(data_path, filename))
                     for filename in dpx_std_file_list]

# These are the wells in the plate that contain the liposome-only solution
# used to calculate the standard curve:
dpx_std_wells = ['A%d' % well_num for well_num in range(1, 13)]

# The list of DPX volumes added at each step
dpx_vol_steps = [0., 1., 2., 4., 6., 8., 10., 10.]

# We know how much of the stock DPX we added at each quenching step, but we
# have to calculate the actual concentrations:
dpx_concs = calc_dpx_concs(dpx_vol_steps)

requench_file_list = [
        '140429_Bax_reDPX_1_0uL.txt',
        '140429_Bax_reDPX_2_1uL.txt',
        '140429_Bax_reDPX_3_2uL.txt',
        '140429_Bax_reDPX_4_4uL.txt',
        '140429_Bax_reDPX_5_6uL.txt',
        '140429_Bax_reDPX_6_8uL.txt',
        '140429_Bax_reDPX_7_10uL.txt'
        ]
requench_file_list = [abspath(join(data_path, filename))
                      for filename in requench_file_list]

# The wells containing the Bid/Bax treated liposomes
bax_requench_wells = [row + col for (col, row) in itertools.product(
                            [str(i) for i in range(1, 13)],
                            #['C'])]
                            ['D', 'E'])]
# The wells containing the Cecropin A treated liposomes
cec_requench_wells = [row + col for (col, row) in itertools.product(
                            [str(i) for i in range(1, 13)],
                            ['C'])]

# The list of DPX volumes added at each step
requench_vol_steps = [0., 1., 2., 4., 6., 8., 10.]

# We know how much of the stock DPX we added at each quenching step, but we
# have to calculate the actual concentrations:
requench_dpx_concs = calc_dpx_concs(requench_vol_steps,
                                    starting_well_vol=100.)

fmax_filename = abspath(join(data_path,
                                '140429_Bax_reDPX_final_triton.txt'))

if __name__ == '__main__':
    plt.ion()

    plt.figure()
    plot_all(bax_wells)

    plt.figure()
    plot_all(cec_wells)

    plot_endpoints_vs_dose(bax_averages, bax_layout)
    plot_endpoints_vs_dose(cec_averages, cec_layout)

    # Calculate the quenching by well
    (i_avgs_by_well, i_sds_by_well, fmax_avg, fmax_sd) = \
            quenching_std_curve_by_well(dpx_std_file_list,
                                        dpx_std_wells, dpx_concs)
    (ka, kd) = fit_std_curve(i_avgs_by_well, i_sds_by_well, dpx_concs)

    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename,
                                         requench_wells, ka, kd,
                                         requench_dpx_concs[-1])

    requenching_analysis(requench_file_list, requench_wells,
                         requench_dpx_concs, fmax_avgs, fmax_sds, ka, kd, 5)

