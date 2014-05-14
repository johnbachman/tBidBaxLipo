from tbidbaxlipo.util.plate_assay import *
from tbidbaxlipo.util.dpx_assay import *
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import tbidbaxlipo.data
from os.path import abspath, join
import collections

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)

timecourse_wells = read_flexstation_kinetics(abspath(join(data_path,
                                             '140425_Bax_DPX_Timecourse.txt')))

# In this experiment, we first prepare liposomes at the same final
# concentration as our wells that will contain Bid and Bax; lyse them with
# Triton; then sequentially add back DPX and measure the quenching at each
# DPX concentration.  Each addition of DPX is read and hence goes in a
# separate file. This list contains the names of files containing the data
# for the standard curve:
dpx_std_file_list = [
    '140425_Bax_DPX_post_triton.txt',
    '140425_Bax_DPX_1_1uL.txt',
    '140425_Bax_DPX_2_1uL.txt',
    '140425_Bax_DPX_3_2uL.txt',
    '140425_Bax_DPX_4_4uL.txt',
    '140425_Bax_DPX_5_4uL.txt',
    '140425_Bax_DPX_6_8uL.txt',
    ]
dpx_std_file_list = [abspath(join(data_path, filename))
                     for filename in dpx_std_file_list]

# These are the wells in the plate that contain the liposome-only solution
# used to calculate the standard curve:
dpx_std_wells = ['A%d' % well_num for well_num in range(1, 13)]

# The list of DPX volumes added at each step
dpx_vol_steps = [0., 1., 1., 2., 4., 4., 8.]

# We know how much of the stock DPX we added at each quenching step, but we
# have to calculate the actual concentrations:
dpx_concs = calc_dpx_concs(dpx_vol_steps)

requench_file_list = [
    '140425_Bax_DPXre_0_0uL.txt',
    '140425_Bax_DPXre_1_1uL.txt',
    '140425_Bax_DPXre_2_2uL.txt',
    '140425_Bax_DPXre_3_4uL.txt',
    '140425_Bax_DPXre_4_6uL.txt',
    '140425_Bax_DPXre_5_8uL.txt']
requench_file_list = [abspath(join(data_path, filename))
                     for filename in requench_file_list]

# The wells containing the Bid/Bax treated liposomes
requench_wells = [row + col for (col, row) in itertools.product(
                            [str(i) for i in range(1, 13)],
                            ['B', 'C', 'D'])]
#requench_wells  = ['B%d' % well_num for well_num in range(1, 13)]
#requench_wells += ['C%d' % well_num for well_num in range(1, 13)]
#requench_wells += ['D%d' % well_num for well_num in range(1, 13)]

# The list of DPX volumes added at each step
requench_vol_steps = [0., 1., 2., 4., 6., 8.]

# We know how much of the stock DPX we added at each quenching step, but we
# have to calculate the actual concentrations:
requench_dpx_concs = calc_dpx_concs(requench_vol_steps,
                                    starting_well_vol=100.)

if __name__ == '__main__':
    plt.ion()

    plt.figure()
    plot_all(timecourse_wells)

    # Calculate the quenching by well
    (i_avgs_by_well, i_sds_by_well, f_max_avg, f_max_sd) = \
                quenching_std_curve_by_well(dpx_std_file_list,
                                            dpx_std_wells, dpx_concs)
    (ka, kd) = fit_std_curve(i_avgs_by_well, i_sds_by_well, dpx_concs)

    requenching_analysis(requench_file_list, requench_wells,
                         requench_dpx_concs, f_max_avg, f_max_sd,
                         ka, kd, 45)
