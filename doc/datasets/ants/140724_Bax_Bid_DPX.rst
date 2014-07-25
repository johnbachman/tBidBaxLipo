.. _140724_Bax_Bid_DPX:

Bax and Bid DPX Requenching Assay (7/24/14)
===========================================

Experimental details
--------------------

.. ipython:: python

    from tbidbaxlipo.plots.layout_140724 import *
    # The first Bax dilution series, Bid 0.625 nM
    print '\n'.join([str(l) for l in bid06_labels])
    # The second Bax dilution series, Bid 2.5 nM
    print '\n'.join([str(l) for l in bid2_labels])
    # The third Bax dilution series, Bid 10 nM
    print '\n'.join([str(l) for l in bid10_labels])
    # The fourth Bax dilution series, Bid 40 nM
    print '\n'.join([str(l) for l in bid40_labels])
    # The wells containing liposomes only:
    print str(dpx_std_wells)
    # The list of DPX volumes added at each std. curve step
    print dpx_vol_steps
    # The list of concentrations at each std. curve step (in Molar)
    print str(dpx_concs)

Raw Timecourses
---------------

.. plot::
    :context:

    from tbidbaxlipo.plots.layout_140724 import *
    plt.close('all')
    plt.figure()
    plot_all(bid06_wells)
    plt.title('Bid 0.625 nM + Bax, timecourses')

    plt.figure()
    plot_all(bid2_wells)
    plt.title('Bid 2.5 nM + Bax, timecourses')

    plt.figure()
    plot_all(bid10_wells)
    plt.title('Bid 10 nM + Bax, timecourses')

    plt.figure()
    plot_all(bid40_wells)
    plt.title('Bid 40 nM + Bax, timecourses')

Endpoints vs. dose
------------------

DPX Background fluorescence
---------------------------

As with 140710, there was an outlier during one of the measurements, which was
discarded in taking the mean.

.. plot::
    :context:

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

Quenching standard curve
------------------------

Here the quenching values are calculated as averages (i.e., the average across
wells at each quenching step, vs. the average of the Triton-only (no DPX)
wells), since there was no Triton in the quenching wells when taking the 0 DPX
measurement.

.. plot::
    :context:

    plt.close('all')
    (i_avgs, i_sds, fmax_avg, fmax_sd) = \
            quenching_std_curve(dpx_std_file_list, dpx_std_wells, dpx_concs,
                                bg_avgs=bg_avgs)
    (ka, kd) = fit_std_curve(i_avgs, i_sds, dpx_concs)
    qd = get_quenching_dict(i_avgs, i_sds, dpx_vols_added)
    final_q = qd[dpx_vols_added[-1]]
    q_outs = np.array(qd.values())

Requenching analysis for Bid 0.625 nM + Bax
-------------------------------------------

.. plot::
    :context:

    plt.close('all')
    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bid06_requench_wells,
                                         final_q, final_bg=bg_avgs[-1])
    requenching_analysis(dpx_std_file_list, bid06_requench_wells, dpx_concs,
                         q_outs, fmax_avgs, fmax_sds, None, None, None, bg_avgs)

Requenching analysis for Bid 2.5 nM + Bax
-----------------------------------------

.. plot::
    :context:

    plt.close('all')
    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bid2_requench_wells,
                                         final_q, final_bg=bg_avgs[-1])
    requenching_analysis(dpx_std_file_list, bid2_requench_wells, dpx_concs,
                         q_outs, fmax_avgs, fmax_sds, None, None, None, bg_avgs)

Requenching analysis for Bid 10 nM + Bax
----------------------------------------

.. plot::
    :context:

    plt.close('all')
    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bid10_requench_wells,
                                         final_q, final_bg=bg_avgs[-1])
    requenching_analysis(dpx_std_file_list, bid10_requench_wells, dpx_concs,
                         q_outs, fmax_avgs, fmax_sds, None, None, None, bg_avgs)

Requenching analysis for Bid 40 nM + Bax
----------------------------------------

.. plot::
    :context:

    plt.close('all')
    (fmax_avgs, fmax_sds) = fmax_by_well(fmax_filename, bid40_requench_wells,
                                         final_q, final_bg=bg_avgs[-1])
    requenching_analysis(dpx_std_file_list, bid40_requench_wells, dpx_concs,
                         q_outs, fmax_avgs, fmax_sds, None, None, None, bg_avgs)

