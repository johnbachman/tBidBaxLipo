from tbidbaxlipo.util.plate_assay import read_flexstation_kinetics, VALUE
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util.error_propagation import calc_ratio_mean_sd
from tbidbaxlipo.util import fitting
import scipy.stats
import itertools

def calc_dpx_concs(dpx_vol_steps, dpx_stock_conc=0.1, starting_well_vol=102.):
    """Calculates the concentration of DPX at each quenching step.

    Parameters
    ----------
    dpx_vol_steps : list of floats
        The amounts (in microliters) of the stock DPX added at each quenching
        step, e.g., [0., 1., 1., 2., 4., 8.] etc.
    dpx_stock_conc : float
        The concentration of the stock DPX solution, in Molar. Defaults to
        0.1 (100 millimolar).
    starting_well_vol : float
        The volume of the well (in microliters) before any DPX is added.
        Defaults to 102 (100 uL well plus 2 uL Triton).

    Returns
    -------
    numpy.array
        The concentrations at each quenching step.
    """

    # The cumulative volumes of stock DPX added at each step
    dpx_vols_added = np.cumsum(dpx_vol_steps)
    # The number of points in the standard curve
    num_dilutions = len(dpx_vols_added)
    dpx_concs = np.zeros(num_dilutions)
    # Calculate the concentrations at each quenching step
    for i in range(num_dilutions):
        total_vol = starting_well_vol + dpx_vols_added[i]
        dpx_concs[i] = (dpx_stock_conc * dpx_vols_added[i]) / total_vol
    return dpx_concs

def quenching_std_curve(dpx_std_file_list, dpx_std_wells, dpx_concs,
                        bg_avgs=None):
    """Calculates and plots the quenching std curve from the raw ANTS data.

    Parameters
    ----------
    dpx_std_file_list : list of strings
        Names of the text files (exported from the FlexStation plate reader)
        containing the raw fluorescence data at each quenching step.
    dpx_std_wells : list of strings
        Names of the replicate wells measured in the quenching assay,
        e.g., ['A1', 'A2', 'A3'] etc.
    dpx_concs : numpy.array
        The concentrations of DPX at each quenching step.

    Returns
    -------
    tuple : (i_avgs, i_sds, fmax_avg)
        A tuple; the first element in the tuple is the array of mean quenching
        values, given as I / I_0, that is, the ratio of the fluorescence at a
        given amount of DPX to the unquenched fluorescence. The second element
        contains the standard deviations of this ratio at each concentration.
        The third element is the average fluorescence intensity of the lysed,
        unquenched liposomes (liposomes after Triton but before DPX addition.
    """
    num_dilutions = len(dpx_concs)

    # We'll store the means and SDs for the fluorescence intensities at
    # each DPX concentration in here:
    dpx_intensity_avgs = np.zeros(num_dilutions)
    dpx_intensity_sds = np.zeros(num_dilutions)

    # Iterate over the files containing the intensity data
    for dilution_index, std_file in enumerate(dpx_std_file_list):
        well_vals = np.array([])
        timecourse_wells = read_flexstation_kinetics(std_file)
        # Iterate over all the wells used in the standard curve, calculate
        # the read average for each one and add it to a list
        for well_name in dpx_std_wells:
            well_vals = np.concatenate((well_vals,
                                        timecourse_wells[well_name][VALUE]))
            #well_avgs.append(np.mean(timecourse_wells[well_name][VALUE]))
        # Now take the average and SD over all of the replicate wells
        dpx_intensity_avgs[dilution_index] = np.mean(well_vals)
        dpx_intensity_sds[dilution_index] = np.std(well_vals)

    # Plot the intensity values at each DPX concentration
    plt.figure()
    plt.errorbar(dpx_concs, dpx_intensity_avgs, yerr=dpx_intensity_sds,
                 color='k', linewidth=2)
    plt.xlabel('[DPX] (M)')
    plt.ylabel('ANTS Fluorescence (RFU)')
    plt.title('ANTS Fluorescence vs. [DPX]')

    # If we have background values, subtract them here
    if bg_avgs is not None:
        if not len(bg_avgs) == num_dilutions:
            return ValueError("The background array must have as many entries "
                              "as there are DPX concentraitons.")
        dpx_intensity_avgs = dpx_intensity_avgs - bg_avgs

    # Now calculate the intensity ratios, I / I_0
    # To get the standard error of the ratio, we have to account for the
    # error in both the numerator and the denominator:
    i_vals = np.zeros(num_dilutions)
    i_sds = np.zeros(num_dilutions)
    for i in range(num_dilutions):
        (i_vals[i], i_sds[i]) = calc_ratio_mean_sd(
                                              dpx_intensity_avgs[i],
                                              dpx_intensity_sds[i],
                                              dpx_intensity_avgs[0],
                                              dpx_intensity_sds[0])
    return (i_vals, i_sds, dpx_intensity_avgs[0], dpx_intensity_sds[0])

def quenching_std_curve_by_well(dpx_std_file_list, dpx_std_wells,
                                dpx_concs):
    """Calculates and plots the quenching std curve from the raw ANTS data.

    Parameters
    ----------
    dpx_std_file_list : list of strings
        Names of the text files (exported from the FlexStation plate reader)
        containing the raw fluorescence data at each quenching step.
    dpx_std_wells : list of strings
        Names of the replicate wells measured in the quenching assay,
        e.g., ['A1', 'A2', 'A3'] etc.
    dpx_concs : numpy.array
        The concentrations of DPX at each quenching step.

    Returns
    -------
    (i_avgs, i_sds)
        A tuple of arrays; the first element in the tuple is the array of mean
        quenching values, given as I / I_0, that is, the ratio of the
        fluorescence at a given amount of DPX to the unquenched fluorescence.
        The second element contains the standard deviations of this ratio at
        each concentration.
    """
    # We know how much of the stock DPX we added at each quenching step, but we
    # have to calculate the actual concentrations:
    num_dilutions = len(dpx_concs)
    num_wells = len(dpx_std_wells)

    # We'll store the means and SDs for the fluorescence intensities at
    # each DPX concentration in here:
    #dpx_intensity_avgs = np.zeros(num_dilutions)
    #dpx_intensity_sds = np.zeros(num_dilutions)
    dpx_well_avgs = np.zeros((num_dilutions, num_wells))
    dpx_well_sds = np.zeros((num_dilutions, num_wells))

    # Iterate over the files containing the intensity data
    for dilution_index, std_file in enumerate(dpx_std_file_list):
        timecourse_wells = read_flexstation_kinetics(std_file)
        # Iterate over all the wells used in the standard curve, calculate
        # the read average for each one and add it to a list
        for well_index, well_name in enumerate(dpx_std_wells):
            dpx_well_avgs[dilution_index, well_index] = \
                            np.mean(timecourse_wells[well_name][VALUE])
            dpx_well_sds[dilution_index, well_index] = \
                            np.std(timecourse_wells[well_name][VALUE])

    # Plot the intensity values at each DPX concentration
    plt.figure()
    for well_index in range(num_wells):
        plt.errorbar(dpx_concs, dpx_well_avgs[:, well_index],
                     yerr=dpx_well_sds[:, well_index])
    plt.xlabel('[DPX] (M)')
    plt.ylabel('ANTS Fluorescence (RFU)')
    plt.title('ANTS Fluorescence vs. [DPX]')

    # Now calculate the intensity ratios, I / I_0, **for each well**
    i_avgs_by_well = np.zeros((num_dilutions, num_wells))
    i_sds_by_well = np.zeros((num_dilutions, num_wells))
    for well_index in range(num_wells):
        for dilution_index in range(num_dilutions):
            # To get the standard error of the ratio, we have to account for
            # the error in both the numerator and the denominator:
            (avg, sd) = calc_ratio_mean_sd(
                                      dpx_well_avgs[dilution_index, well_index],
                                      dpx_well_sds[dilution_index, well_index],
                                      dpx_well_avgs[0, well_index],
                                      dpx_well_sds[0, well_index])
            i_avgs_by_well[dilution_index, well_index] = avg
            i_sds_by_well[dilution_index, well_index] = sd

    # Plot the individual I / I_0 curves
    plt.figure()
    for well_index in range(num_wells):
        plt.errorbar(dpx_concs, i_avgs_by_well[:, well_index],
                     yerr=i_sds_by_well[:, well_index])
    plt.xlabel('[DPX] (M)')
    plt.ylabel('$I / I_0$')
    plt.title('Quenching vs. [DPX]')

    # Now calculate the I / I_0 curve averaged over all replicates NOTE: we're
    # not taking into account the individual variances at each well (which are
    # relatively small) when we are calculate these means/sds. Doing so
    # could represent a slight improvement.
    i_avgs = np.mean(i_avgs_by_well, axis=1)
    i_sds = np.std(i_avgs_by_well, axis=1)

    fmax_avg = np.mean(dpx_well_avgs[0, :])
    fmax_sd = np.std(dpx_well_avgs[0, :])

    return (i_avgs, i_sds, fmax_avg, fmax_sd)

def fit_std_curve(i_vals, i_sds, dpx_concs):
    # Plot the intensity ratios as a function of [DPX]:
    plt.figure()
    plt.errorbar(dpx_concs, i_vals, yerr=i_sds, linestyle='', marker='o',
                 color='k', linewidth=2)
    plt.xlabel('[DPX] (M)')
    plt.ylabel('$I / I_0$')
    plt.title('Quenching vs. [DPX]')

    # Fit the curve using scipy least squares
    kd = fitting.Parameter(0.05)
    ka = fitting.Parameter(0.490)
    def eq2(dpx):
        return 1. / ((1 + kd()*dpx) * (1 + ka()*dpx))
    fitting.fit(eq2, [kd, ka], i_vals, dpx_concs)
    # Plot curve fit
    dpx_interp = np.linspace(dpx_concs[0], dpx_concs[-1], 100)
    plt.plot(dpx_interp, eq2(dpx_interp), color='r', linewidth=2)
    plt.xlabel('[DPX] (mM)')
    plt.ylabel('$I / I_0$')
    plt.title('DPX quenching, standard curve')
    print "kd: %f" % kd()
    print "ka: %f" % ka()

    return (ka(), kd())

def quenching_func(ka, kd, dpx_concs):
    return (1 + kd * dpx_concs) * (1 + ka * dpx_concs)

def fit_std_curve_by_pymc(i_vals, i_sds, dpx_concs):
    import pymc
    from pymc import Uniform, stochastic, deterministic, MCMC
    from pymc import Matplot
    # Define prior distributions for both Ka and Kd
    ka = Uniform('ka', lower=0, upper=1000)
    kd = Uniform('kd', lower=0, upper=1000)

    @stochastic(plot=True, observed=True)
    def quenching_model(ka=ka, kd=kd, value=i_vals):
        pred_i = quenching_func(ka, kd, dpx_concs)
        # The first concentration in dpx_concs should always be zero
        # (that is, the first point in the titration should be the
        # unquenched fluorescence), so we assert that here:
        assert dpx_concs[0] == 0
        # The reason this is necessary is that in the likelihood calculation
        # we skip the error for the first point, since (when the std. err
        # is calculated by well) the error is 0 (the I / I_0 ratio is
        # always 1 for each well, the the variance/SD across the wells is 0).
        # If we don't skip this first point, we get nan for the likelihood.
        # In addition, the model always predicts 1 for the I / I_0 ratio
        # when the DPX concentration is 0, so it contributes nothing to
        # the overall fit.
        return -np.sum((value[1:] - pred_i[1:])**2 / (2 * i_sds[1:]**2))

    pymc_model = pymc.Model([ka, kd, quenching_model])
    mcmc = MCMC(pymc_model)
    mcmc.sample(iter=155000, burn=5000, thin=150)
    Matplot.plot(mcmc)

    plt.figure()
    num_to_plot = 1000
    ka_vals = mcmc.trace('ka')[:]
    kd_vals = mcmc.trace('kd')[:]
    if num_to_plot > len(ka_vals):
        num_to_plot = len(ka_vals)
    for i in range(num_to_plot):
        plt.plot(dpx_concs, quenching_func(ka_vals[i], kd_vals[i], dpx_concs),
                 alpha=0.01, color='r')
    plt.errorbar(dpx_concs, i_vals, yerr=i_sds, linestyle='', marker='o',
            color='k', linewidth=2)

    return (ka_vals, kd_vals)

def requenching_analysis(requench_file_list, requench_wells,
                         requench_dpx_concs, q_outs, fmax_avg, fmax_sd,
                         ka, kd, dpx_0, bg_avgs=None, do_plot=True):
    """Calculates and plots the quenching std curve from the raw ANTS data.

    Parameters
    ----------
    requench_file_list : list of strings
        Names of the text files (exported from the FlexStation plate reader)
        containing the raw fluorescence data at each quenching step.
    requench_wells : list of strings
        Names of the replicate wells measured in the quenching assay,
        e.g., ['A1', 'A2', 'A3'] etc.
    requench_dpx_concs : numpy.array
        The concentrations of DPX at each quenching step.
    """
    # FIXME add other comments to docstring

    # We know how much of the stock DPX we added at each quenching step, but we
    # have to calculate the actual concentrations:
    num_dilutions = len(requench_dpx_concs)
    num_wells = len(requench_wells)

    # We'll store the means and SDs for the fluorescence intensities at
    # each DPX concentration in here:
    requench_well_avgs = np.zeros((num_dilutions, num_wells))
    requench_well_sds = np.zeros((num_dilutions, num_wells))

    # Iterate over the files containing the intensity data
    for dilution_index, requench_file in enumerate(requench_file_list):
        timecourse_wells = read_flexstation_kinetics(requench_file)
        # Iterate over all the wells, calculate the read average for each one
        # and add it to a list
        for well_index, well_name in enumerate(requench_wells):
            if bg_avgs is not None:
                well_values = (timecourse_wells[well_name][VALUE] -
                               bg_avgs[dilution_index])
            else:
                well_values = timecourse_wells[well_name][VALUE]

            requench_well_avgs[dilution_index, well_index] = \
                            np.mean(well_values)
            requench_well_sds[dilution_index, well_index] = \
                            np.std(well_values)

    # Plot the intensity values at each DPX concentration
    if do_plot:
        plt.figure()
        for well_index in range(num_wells):
            plt.errorbar(requench_dpx_concs, requench_well_avgs[:, well_index],
                         yerr=requench_well_sds[:, well_index])
        plt.xlabel('[DPX] (M)')
        plt.ylabel('ANTS Fluorescence (RFU)')
        plt.title('ANTS Fluorescence vs. [DPX]')

    # Now calculate the total quenching, I / F_max, for each well
    q_tot_by_well_avgs = np.zeros((num_dilutions, num_wells))
    q_tot_by_well_sds = np.zeros((num_dilutions, num_wells))
    for well_index, well_name in enumerate(requench_wells):
        # If we've been given a list of fmax avgs/sds for each well, use the
        # one for this well; otherwise, assume we've been given a single value
        # for all wells and use that
        try:
            cur_fmax_avg = fmax_avg[well_index]
            cur_fmax_sd = fmax_sd[well_index]
        except IndexError:
            cur_fmax_avg = fmax_avg
            cur_fmax_sd = fmax_sd

        # The error associated with the total quenching depends on the
        # error of the reads for this well/dilution, and the error of the
        # maximum fluorescence
        for dilution_index in range(num_dilutions):
            (q_tot_avg, q_tot_sd) = calc_ratio_mean_sd(
                        requench_well_avgs[dilution_index, well_index],
                        requench_well_sds[dilution_index, well_index],
                        cur_fmax_avg, cur_fmax_sd)
            q_tot_by_well_avgs[dilution_index, well_index] = q_tot_avg
            q_tot_by_well_sds[dilution_index, well_index] = q_tot_sd

    # Calculate the (predicted) Q_out for each of the DPX concentrations used
    #q_outs = np.array([1. / quenching_func(ka, kd, dpx_conc)
    #          for dpx_conc in requench_dpx_concs])
    # FIXME have q_out get associated error by sampling from ka/kd dists

    # Plot the Q_total / Q_out curves (should be straight lines)
    q_ins = np.zeros(num_wells)
    q_in_errs = np.zeros(num_wells)
    f_outs = np.zeros(num_wells)
    f_out_errs = np.zeros(num_wells)
    #plt.figure()
    for well_index, well_name in enumerate(requench_wells):
        #if well_index < 30:
        #    continue
        if do_plot:
            plt.figure()
            #plt.subplot(6, 6, well_index+1)
            plt.errorbar(q_outs, q_tot_by_well_avgs[:, well_index],
                         q_tot_by_well_sds[:, well_index],
                         linestyle='', color='k', linewidth=2)
            plt.title('$Q_{total}$ vs. $Q_{out}$ for well %s' % well_name)
            plt.xlim((0, 1))
            plt.ylim((0, 1))
            plt.xlabel('$Q_{out}$')
            plt.ylabel('$Q_{total}$')
            # Add tickmarks only to the plots on the edges
            #if well_index + 1 > 30 
            #    plt.xticks([0, 0.5, 1.0])
            #    plt.xlabel('$Q_{out}$')
            #else:
            #    plt.xticks([])
            #if well_index % 6 == 0:
            #    plt.yticks([0, 0.5, 1.0])
            #    plt.ylabel('$Q_{total}$')
            #else:
            #    plt.yticks([])

        linfit = scipy.stats.linregress(q_outs,
                                        q_tot_by_well_avgs[:, well_index])
        f_out = linfit[0]
        intercept = linfit[1]
        #std_err = linfit[4]
        n = float(len(q_outs))
        sd_slope = linfit[4]
        sd_intercept = np.sqrt((sd_slope ** 2) * np.sum(q_outs ** 2) / n)

        # Calculate the internal quenching
        (q_in_avg, q_in_sd) = calc_ratio_mean_sd(intercept, sd_intercept,
                                                 1 - f_out, sd_slope)

        # THIS IS ALL WRONG!
        #mean_x = np.mean(q_outs)
        #sx2 =  np.sum((q_outs - mean_x)**2)
        #sd_intercept = std_err * np.sqrt(1. / len(q_outs) + mean_x*mean_x/sx2)
        #sd_slope = std_err * np.sqrt(1./sx2)
        if do_plot:
            plt.plot(q_outs, (q_outs * f_out) + intercept, color='r',
                     linewidth=1)
            plt.text(0.05, 0.9, '$F_{out}\ \mathrm{(slope)}\ = %f$' % f_out)
            plt.text(0.05, 0.8, '$Q_{in}(1 - F_{out})\ \mathrm{(intercept)}\ = %f$'
                                % intercept)
            plt.text(0.05, 0.7, '$Q_{in} = %f$' % q_in_avg)
            plt.text(0.05, 0.6, '$R^2 = %f$' % linfit[2])

        # Save the results for this well
        f_outs[well_index] = f_out
        f_out_errs[well_index] = sd_slope
        q_ins[well_index] = q_in_avg
        q_in_errs[well_index] = q_in_sd
    #plt.tight_layout()

    # This error is not meaningful because it doesn't account for errorbars!
    # FIXME FIXME FIXME
    linfit = scipy.stats.linregress(f_outs, q_ins)
    f_out_pred = np.linspace(0, 1, 100)

    if do_plot:
        plt.figure()
        plt.errorbar(f_outs, q_ins, xerr=f_out_errs, yerr=q_in_errs,
                     marker='o', color='k', linestyle='')
        #plt.plot(f_outs[0:12], q_ins[0:12], marker='o', color='r',
        #         linestyle='')
        #plt.plot(f_outs[12:24], q_ins[12:24], marker='o', color='g',
        #         linestyle='')
        #plt.plot(f_outs[24:36], q_ins[24:36], marker='o', color='b',
        #         linestyle='')
        plt.xlabel('$F_{out}$')
        plt.ylabel('$Q_{in}$')
        plt.xlim([-0.2, 1.1])
        plt.ylim([-0.2, 1.1])
        plt.title('$Q_{in}$ vs. $F_{out}$')

        plt.plot(f_out_pred, (f_out_pred * linfit[0]) + linfit[1], color='r')
        plt.text(0.1, 0.9, 'Slope = %f' % linfit[0])
        plt.text(0.1, 0.85, 'Intercept = %f' % linfit[1])
        plt.text(0.1, 0.8, '$R^2 = %f$' % linfit[2])
        plt.text(0.1, 0.75, 'p = %f' % linfit[3])

    return {'f_outs': f_outs,
            'q_ins': q_ins,
            'f_out_errs': f_out_errs,
            'q_in_errs': q_in_errs}


    #alpha = fitting.Parameter(1.)
    #dpx_0 = fitting.Parameter(18.1)
    #dpx_0 = fitting.Parameter(0.0144)
    #def q_in_func(f_out):
    #    return (1. /
    #            (((1 + kd * dpx_0()) * (1 - f_out) ** alpha()) *
    #             ((1 + ka * dpx_0()) * (1 - f_out) ** alpha())))
    #fitting.fit(q_in_func, [alpha, dpx_0], q_ins, f_outs)
    #fitting.fit(q_in_func, [alpha, dpx_0], q_ins, f_outs)
    #plt.plot(f_out_pred, q_in_func(f_out_pred), linewidth=2, color='r')
    #print "dpx_0: %f" % dpx_0()
    #print "alpha: %f" % alpha()

    #for (dpx_val, alpha_val) in itertools.product((0.005, 0.002),
    #                                              (alpha(), 0)):
    #    alpha.set(alpha_val)
    #    dpx_0.set(dpx_val)
    #    plt.plot(f_out_pred, q_in_func(f_out_pred), linewidth=2, color='b')

def fmax_by_well(fmax_file, requench_wells, final_q, final_bg=None):
    num_wells = len(requench_wells)

    # We'll store the means and SDs for the fluorescence intensities here
    fmax_avgs = np.zeros(num_wells)
    fmax_sds = np.zeros(num_wells)

    timecourse_wells = read_flexstation_kinetics(fmax_file)
    # Iterate over all the wells, calculate the read average for each one
    # and add it to a list
    for well_index, well_name in enumerate(requench_wells):
        quenched_vals = timecourse_wells[well_name][VALUE]
        # First, if possible, subtract the background from the raw value
        if final_bg is not None:
            quenched_vals -= final_bg
        # Now we calculate what the values would be if they were unquenched
        fmax_vals = quenched_vals / final_q
        fmax_avgs[well_index] = np.mean(fmax_vals)
        fmax_sds[well_index] = np.std(fmax_vals)

    return (fmax_avgs, fmax_sds)

def get_quenching_dict(i_avgs, i_sds, dpx_vols_added):
    return dict([(dpx_vols_added[i], i_avgs[i])
                 for i in range(len(dpx_vols_added))])

