from tbidbaxlipo.util.dpx_assay import *
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import tbidbaxlipo.data

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)

if __name__ == '__main__':
    plt.ion()

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
    dpx_std_file_list = [os.path.abspath(os.path.join(data_path, filename))
                         for filename in dpx_std_file_list]

    # These are the wells in the plate that contain the liposome-only solution
    # used to calculate the standard curve:
    dpx_std_wells = ['A%d' % well_num for well_num in range(1, 13)]

    # The list of DPX volumes added at each step
    dpx_vol_steps = [0., 1., 2., 4., 6., 8., 10., 10.]

    # We know how much of the stock DPX we added at each quenching step, but we
    # have to calculate the actual concentrations:
    dpx_concs = calc_dpx_concs(dpx_vol_steps)

    # Calculate the quenching by well
    (i_avgs_by_well, i_sds_by_well, fmax_avg, fmax_sd) = \
            quenching_std_curve_by_well(dpx_std_file_list,
                                        dpx_std_wells, dpx_concs)
    (ka, kd) = \
            fit_std_curve(i_avgs_by_well, i_sds_by_well, dpx_concs)

    """
    (ka_vals_by_well, kd_vals_by_well) = \
                fit_std_curve_by_pymc(i_avgs_by_well, i_sds_by_well, dpx_concs)
    print "Mean for posterior K_D: %f" % np.mean(kd_vals_by_well)
    print "SD for posterior K_D: %f" % np.std(kd_vals_by_well)
    print "Mean for posterior K_A: %f" % np.mean(ka_vals_by_well)
    print "SD for posterior K_A: %f" % np.std(ka_vals_by_well)

    # Now do it where the quenching is not calculated by well
    (i_avgs, i_sds, fmax_avg, fmax_sd) = quenching_std_curve(dpx_std_file_list,
                                           dpx_std_wells, dpx_concs)
    (ka, kd) = fit_std_curve(i_avgs, i_sds, dpx_concs)

    #(ka_vals, kd_vals) = fit_std_curve_by_pymc(i_avgs, i_sds, dpx_concs)
    #print "Mean for posterior K_D: %f" % np.mean(kd_vals)
    #print "SD for posterior K_D: %f" % np.std(kd_vals)
    #print "Mean for posterior K_A: %f" % np.mean(ka_vals)
    #print "SD for posterior K_A: %f" % np.std(ka_vals)
    """
    requench_file_list = [
            '140429_Bax_reDPX_1_0uL.txt',
            '140429_Bax_reDPX_2_1uL.txt',
            '140429_Bax_reDPX_3_2uL.txt',
            '140429_Bax_reDPX_4_4uL.txt',
            '140429_Bax_reDPX_5_6uL.txt',
            '140429_Bax_reDPX_6_8uL.txt',
            '140429_Bax_reDPX_7_10uL.txt'
            ]
    requench_file_list = [os.path.abspath(os.path.join(data_path, filename))
                          for filename in requench_file_list]

    # The wells containing the Bid/Bax treated liposomes
    requench_wells = [row + col for (col, row) in itertools.product(
                                [str(i) for i in range(1, 13)],
                                #['C'])]
                                ['D', 'E'])]

    # The list of DPX volumes added at each step
    requench_vol_steps = [0., 1., 2., 4., 6., 8., 10.]

    # We know how much of the stock DPX we added at each quenching step, but we
    # have to calculate the actual concentrations:
    requench_dpx_concs = calc_dpx_concs(requench_vol_steps,
                                        starting_well_vol=100.)

    fmax_filename = os.path.abspath(os.path.join(data_path,
                                    '140429_Bax_reDPX_final_triton.txt'))
    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename,
                                         requench_wells, ka, kd,
                                         requench_dpx_concs[-1])

    requenching_analysis(requench_file_list, requench_wells,
                         requench_dpx_concs, fmax_avgs, fmax_sds, ka, kd, 5)
