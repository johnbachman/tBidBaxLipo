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

# Calculate starting dose for first Bax dilution series
bax_stock_conc = 5900. # nM
total_vol = 150. # uL
bax_stock1 = 50. # uL
bax_stock2 = 41.5 # uL
bax_start_conc1 = (bax_stock_conc * bax_stock1)/total_vol
bax_start_conc2 = (bax_stock_conc * bax_stock2)/total_vol

bax_labels1 = dose_series_replicate_list('Bax', bax_start_conc1, 2/3.,
                                         num_doses=12,
                                         start_row='B', end_row='B',
                                         start_col=1, end_col=12)
bax_labels2 = dose_series_replicate_list('Bax', bax_start_conc2, 2/3.,
                                         num_doses=12,
                                         start_row='C', end_row='C',
                                         start_col=1, end_col=12)
bax_layout = collections.OrderedDict(bax_labels1 + bax_labels2)

cec_layout = collections.OrderedDict(dose_series_replicate_list(
                    'CecropinA', 20000, 2/3., num_doses=12,
                    start_row='D', end_row='D', start_col=1, end_col=12))

timecourse_file = abspath(join(data_path, '140514_Bax_DPX_43C_timecourse.txt'))
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
        '140514_Bax_DPX_43C_1_0uL.txt',
        '140514_Bax_DPX_43C_2_1uL.txt',
        '140514_Bax_DPX_43C_3_2uL.txt',
        '140514_Bax_DPX_43C_4_4uL.txt',
        '140514_Bax_DPX_43C_5_6uL.txt',
        '140514_Bax_DPX_43C_6_8uL.txt',
        '140514_Bax_DPX_43C_7_10uL.txt']
dpx_std_file_list = [abspath(join(data_path, filename))
                     for filename in dpx_std_file_list]

# These are the wells in the plate that contain the liposome-only solution
# used to calculate the standard curve:
dpx_std_wells = ['A%d' % well_num for well_num in range(1, 13)]

# The list of DPX volumes added at each step
dpx_vol_steps = [0., 1., 2., 4., 6., 8., 10.]

# The cumulative volumes of stock DPX added at each step
dpx_vols_added = np.cumsum(dpx_vol_steps)

# We know how much of the stock DPX we added at each quenching step, but we
# have to calculate the actual concentrations:
dpx_concs = calc_dpx_concs(dpx_vol_steps)

requench_file_list = [
        '140514_Bax_DPX_43C_1_0uL.txt',
        '140514_Bax_DPX_43C_2_1uL.txt',
        '140514_Bax_DPX_43C_3_2uL.txt',
        '140514_Bax_DPX_43C_4_4uL.txt',
        '140514_Bax_DPX_43C_5_6uL.txt',
        '140514_Bax_DPX_43C_6_8uL.txt',
        '140514_Bax_DPX_43C_7_10uL.txt']
requench_file_list = [abspath(join(data_path, filename))
                      for filename in requench_file_list]

# The wells containing the Bid/Bax treated liposomes
bax_requench_wells = [row + col for (col, row) in itertools.product(
                            [str(i) for i in range(1, 13)],
                            #['C'])]
                            ['B', 'C'])]
# The wells containing the Cecropin A treated liposomes
cec_requench_wells = [row + col for (col, row) in itertools.product(
                            [str(i) for i in range(1, 13)],
                            ['D'])]

# The list of DPX volumes added at each step
requench_vol_steps = [0., 1., 2., 4., 6., 8., 10.]

# We know how much of the stock DPX we added at each quenching step, but we
# have to calculate the actual concentrations:
requench_dpx_concs = calc_dpx_concs(requench_vol_steps,
                                    starting_well_vol=100.)

fmax_filename = abspath(join(data_path, '140514_Bax_DPX_43C_triton.txt'))

if __name__ == '__main__':
    plt.ion()

    #plt.figure()
    #plot_all(bax_wells)

    #plt.figure()
    #plot_all(cec_wells)

    #plot_endpoints_vs_dose(bax_averages, bax_layout)
    #plot_endpoints_vs_dose(cec_averages, cec_layout)

    # Calculate the quenching by well
    (i_avgs_by_well, i_sds_by_well, fmax_avg, fmax_sd) = \
            quenching_std_curve_by_well(dpx_std_file_list,
                                        dpx_std_wells, dpx_concs)
    qd = get_quenching_dict(i_avgs_by_well, i_sds_by_well, dpx_vols_added)

    final_q = qd[dpx_vols_added[-1]]

    #(ka, kd) = fit_std_curve(i_avgs_by_well, i_sds_by_well, dpx_concs)

    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bax_requench_wells,
                                         final_q)
    q_outs = np.array(qd.values())

    requenching_analysis(requench_file_list, bax_requench_wells,
                         requench_dpx_concs, q_outs, fmax_avgs, fmax_sds,
                         None, None, None)

