from tbidbaxlipo.util.dpx_assay import *
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import tbidbaxlipo.data
from tbidbaxlipo.util.plate_assay import *
from os.path import abspath, join
import collections
from tbidbaxlipo.plots import titration_fits as tf
import matplotlib.patches as patches
from tbidbaxlipo.util import fitting
from scipy import stats

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

    # 0.6 nM Bid
    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bid06_requench_wells,
                                         final_q, final_bg=bg_avgs[-1])
    req_results = \
          requenching_analysis(dpx_std_file_list, bid06_requench_wells, dpx_concs,
            q_outs, fmax_avgs, fmax_sds, None, None, None, bg_avgs,
            do_plot=False)
    plot_requenching_result(req_results, '140724_requench_bid_06')

    # 2 nM Bid
    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bid2_requench_wells,
                                         final_q, final_bg=bg_avgs[-1])
    req_results = \
          requenching_analysis(dpx_std_file_list, bid2_requench_wells, dpx_concs,
            q_outs, fmax_avgs, fmax_sds, None, None, None, bg_avgs,
            do_plot=False)
    plot_requenching_result(req_results, '140724_requench_bid_2')

    # 10 nM cBid
    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bid10_requench_wells,
                                         final_q, final_bg=bg_avgs[-1])
    req_results = \
          requenching_analysis(dpx_std_file_list, bid06_requench_wells, dpx_concs,
            q_outs, fmax_avgs, fmax_sds, None, None, None, bg_avgs,
            do_plot=False)
    plot_requenching_result(req_results, '140724_requench_bid_10')

    # 40 nM cBid
    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bid40_requench_wells,
                                         final_q, final_bg=bg_avgs[-1])
    req_results = \
          requenching_analysis(dpx_std_file_list, bid40_requench_wells, dpx_concs,
                       q_outs, fmax_avgs, fmax_sds, None, None, None, bg_avgs,
                       do_plot=False)
    plot_requenching_result(req_results, '140724_requench_bid_40')


def bid_bax_kinetic_analysis(bid_wells, bid_layout, bid_conc_str):
    # Improve this by calculating the Fmax value based on Triton/DPX

    # First normalize by the background
    # Background is no Bax condition.
    bg_key_str = 'Bid %s nM, Bax 0.0 nM' % bid_conc_str
    bg_tc = bid_wells[bid_layout[bg_key_str][0]]
    bg_time = bg_tc[TIME]
    bg_value = bg_tc[VALUE]

    fmax_arr = []
    k1_arr = []
    k2_arr = []
    bax_concs = []
    plt.figure(1)
    for conc_name, well_list in bid_layout.iteritems():
        bax_conc = conc_name.split(' ')[4]
        bax_concs.append(bax_conc)
        if bax_conc == 0:
            continue
        # For this expt, there is only one well per condition (no replicates)
        well_name = well_list[0]
        time = bid_wells[well_name][TIME] + 240
        value = bid_wells[well_name][VALUE]
        bg_corr_value = value / bg_value

        # In lieu of triton, assume full perm
        maxval = 113.
        max_norm_val = maxval / bg_value[-1]
        norm_value = (bg_corr_value - 1.0) / (max_norm_val - 1.0)

        fit = tf.OneExpFmax()
        karr = fit.fit_timecourse(time, norm_value)
        k1_arr.append(karr[0])
        fmax_arr.append(karr[1])
        #k2_arr.append(k2)
        plt.plot(time, norm_value)
        plt.plot(time, fit.fit_func(time, karr), color='r')

    plt.figure('Titration')
    plt.subplot(1,2,1)
    plt.plot(bax_concs, fmax_arr, marker='o')
    plt.title('Fmax')
    plt.subplot(1,2,2)
    plt.plot(bax_concs, k1_arr, marker='o')
    plt.title('k1')
    #plt.figure()
    #plt.plot(bax_concs, k2_arr, marker='o')
    #plt.title('k2')


    # Calculate pore sizes
    set_fig_params_for_publication()
    plt.ion()
    #(fmax_arr, k1_arr, k2_arr, conc_list) = \
    #     titration_fits.plot_two_exp_fits(bgsub_norm_wells, layout, plot=False)

    conc_lipos = 2.17
    conc_list = np.array([float(c) for c in bax_concs])
    fmax_means = fmax_arr
    #fmax_stds = np.std(fmax_arr, axis=0)
    ratios = conc_list / conc_lipos
    isf_pore_sizes = []
    k_max = 20
    for i, ratio in enumerate(ratios):
        if i == 0:
            continue
        # Probability of permeabilization is
        # 1 - poisson.cdf(pore_size, ratio)
        pore_sizes = np.arange(1, k_max)
        probs = [1 - stats.poisson.cdf(pore_size - 1, ratio)
                 for pore_size in pore_sizes]
        isf_pore_size = stats.poisson.isf(fmax_means[i], ratio) + 1
        isf_pore_sizes.append(isf_pore_size)

        """
        plt.figure()
        plt.plot(pore_sizes, probs, marker='o')
        plt.xlabel('Pore size')
        plt.ylim([-0.05, 1.05])
        plt.title('Ratio of %f, isf pore size: %f' % (ratio, isf_pore_size))
        sd_lines = [fmax_means[i] + fmax_stds[i]*num_sds
                    for num_sds in np.arange(-2, 3)]
        plt.hlines(sd_lines, 0, k_max, color='r')
        """

    # Plot min pore size vs. Bax, with linear fit
    fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
    ax = fig.gca()
    ax.plot(ratios[1:], isf_pore_sizes, marker='o', markersize=2,
            linestyle='') # Skip the 0 Bax pt
    ax.set_xlabel('[Bax]/[Lipo]')
    ax.set_ylabel('Predicted pore size')
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    lbound = 0.1
    ubound = 200
    #ax.set_xlim(lbound, 1000)
    #ax.set_ylim(1, 200)
    ratios=ratios[1:]
    lin_fit = stats.linregress(ratios[:3], isf_pore_sizes[:3])
    slope = lin_fit[0]
    intercept = lin_fit[1]
    print intercept
    ax.plot(ratios[1:], slope*ratios[1:] + intercept, color='r')
    interp_range = np.linspace(lbound, ratios[1])
    ax.plot(interp_range, slope*interp_range + intercept, color='r',
            linestyle='--', dashes=(2, 2))
    ax.plot(interp_range, [slope*ratios[1] + intercept] * len(interp_range),
            color='gray', linestyle='--', dashes=(2, 2))
    #ax.set_title('Slope %f, intercept %f' % (slope, intercept),
    #             fontsize=fontsize)
    fig.subplots_adjust(left=0.22, bottom=0.19)
    format_axis(ax)


def bid_bax_pore_analysis(bid_wells, bid_layout, bid_conc_str):
    # First normalize by the background
    # Background is no Bax condition.
    bg_key_str = 'Bid %s nM, Bax 0.0 nM' % bid_conc_str
    bg_tc = bid_wells[bid_layout[bg_key_str][0]]
    bg_time = bg_tc[TIME]
    bg_value = bg_tc[VALUE]

    fmax_arr = []
    k1_arr = []
    k2_arr = []
    bax_concs = []
    plt.figure()
    for conc_name, well_list in bid_layout.iteritems():
        bax_conc = conc_name.split(' ')[4]
        bax_concs.append(bax_conc)
        # For this expt, there is only one well per condition (no replicates)
        well_name = well_list[0]
        time = bid_wells[well_name][TIME] + 240
        value = bid_wells[well_name][VALUE]
        bg_corr_value = value / bg_value

        # In lieu of triton, assume full perm
        maxval = 113.
        max_norm_val = maxval / bg_value[-1]
        norm_value = (bg_corr_value - 1.0) / (max_norm_val - 1.0)

        # Now calculate pores
        pores = -np.log(1 - norm_value)
        plt.plot(time, pores)

        fit = tf.LinkedEq(bax_conc, initial_guesses=[1e-2, 1.25e-4, 2])
        karr = fit.fit_timecourse(time, pores)
        k1_arr.append(karr[0])
        k2_arr.append(karr[1])
        fmax_arr.append(karr[2])
        plt.plot(time, fit.fit_func(time, karr), color='r')

    plt.figure('Titration')
    plt.subplot(1,3,1)
    plt.plot(bax_concs, fmax_arr, marker='o')
    plt.title('Fmax')
    plt.subplot(1,3,2)
    plt.plot(bax_concs, k1_arr, marker='o')
    plt.title('k1')
    plt.subplot(1,3,3)
    plt.plot(bax_concs, k2_arr, marker='o')
    plt.title('k2')
    #plt.figure()
    #plt.plot(bax_concs, k2_arr, marker='o')
    #plt.title('k2')

if __name__ == '__main__':
    main()

    #bid_bax_pore_analysis(bid06_wells, bid06_layout, '0.63')
    #bid_bax_pore_analysis(bid40_wells, bid40_layout, '40')

    #bid_bax_kinetic_analysis(bid06_wells, bid06_layout, '0.63')
    #bid_bax_kinetic_analysis(bid2_wells, bid2_layout, '2.5')
    #bid_bax_kinetic_analysis(bid10_wells, bid10_layout, '10')
    #bid_bax_kinetic_analysis(bid40_wells, bid40_layout, '40')
