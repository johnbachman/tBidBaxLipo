from tbidbaxlipo.util.plate_assay import *
import itertools
import pickle
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import fitting, set_fig_params_for_publication, emcee_fit
from tbidbaxlipo.plots.titration_fits import TwoExp, OneExpFmax
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from tbidbaxlipo.models.nbd import multiconf
from tbidbaxlipo.models import one_cpt

layout = collections.OrderedDict([
        ('Bax 185 nM, Lipos 1 mg/ml',  ['A1']),
        ('Bax 185 nM, Lipos 0.5 mg/ml',  ['A2']),
        ('Bax 185 nM, Lipos 0.25 mg/ml',  ['A3']),
        ('Bax 185 nM, Lipos 0.125 mg/ml',  ['A4']),
        ('Bax 185 nM, Lipos 0.063 mg/ml',  ['A5']),
        ('Bax 185 nM, Lipos 0.031 mg/ml',  ['A6']),
        ('Bax 185 nM, Lipos 0.016 mg/ml',  ['A7']),
        ('Bax 185 nM, Lipos 0.008 mg/ml',  ['A8']),
        ('Bax 185 nM, Lipos 0.004 mg/ml',  ['A9']),
        ('Bax 185 nM, Lipos 0.002 mg/ml',  ['A10']),
        ('Bax 185 nM, Lipos 0.001 mg/ml',  ['A11']),
        ('Bax 185 nM, Lipos 0 mg/ml',  ['A12']),
        ('Lipos 1 mg/ml',  ['B1']),
        ('Lipos 0.5 mg/ml',  ['B2']),
        ('Lipos 0.25 mg/ml',  ['B3']),
        ('Lipos 0.125 mg/ml',  ['B4']),
        ('Lipos 0.063 mg/ml',  ['B5']),
        ('Lipos 0.031 mg/ml',  ['B6']),
        ('Lipos 0.016 mg/ml',  ['B7']),
        ('Lipos 0.008 mg/ml',  ['B8']),
        ('Lipos 0.004 mg/ml',  ['B9']),
        ('Lipos 0.002 mg/ml',  ['B10']),
        ('Lipos 0.001 mg/ml',  ['B11']),
        ('Lipos 0 mg/ml',  ['B12']),
    ])
data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '140318_NBD_Bax_BimBH3_lipo_titration.txt'))

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_flexstation_kinetics(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = extract(wells_to_read, timecourse_wells)

lipo_bg_conditions = [
        'Lipos 1 mg/ml',
        'Lipos 0.5 mg/ml',
        'Lipos 0.25 mg/ml',
        'Lipos 0.125 mg/ml',
        'Lipos 0.063 mg/ml',
        'Lipos 0.031 mg/ml',
        'Lipos 0.016 mg/ml',
        'Lipos 0.008 mg/ml',
        'Lipos 0.004 mg/ml',
        'Lipos 0.002 mg/ml',
        'Lipos 0.001 mg/ml',
        'Lipos 0 mg/ml']
lipo_bg_layout = extract(lipo_bg_conditions, layout)
lipo_bg_wells = extract([layout[cond][0] for cond in lipo_bg_conditions],
                        timecourse_wells)

bax_lipo_conditions = [
        'Bax 185 nM, Lipos 1 mg/ml',
        'Bax 185 nM, Lipos 0.5 mg/ml',
        'Bax 185 nM, Lipos 0.25 mg/ml',
        'Bax 185 nM, Lipos 0.125 mg/ml',
        'Bax 185 nM, Lipos 0.063 mg/ml',
        'Bax 185 nM, Lipos 0.031 mg/ml',
        'Bax 185 nM, Lipos 0.016 mg/ml',
        'Bax 185 nM, Lipos 0.008 mg/ml',
        'Bax 185 nM, Lipos 0.004 mg/ml',
        'Bax 185 nM, Lipos 0.002 mg/ml',
        'Bax 185 nM, Lipos 0.001 mg/ml',
        'Bax 185 nM, Lipos 0 mg/ml', ]
bax_lipo_layout = extract(bax_lipo_conditions, layout)
bax_lipo_wells = extract([layout[cond][0] for cond in bax_lipo_conditions],
                         timecourse_wells)

# Normalized and background subtracted
bgsub_wells = subtract_background_set(bax_lipo_wells, lipo_bg_wells)

(bgsub_averages, bgsub_sds) = averages(bgsub_wells, bax_lipo_layout)

#bgsub_norm_wells = subtract_background(norm_wells, background)

# Normalized, background subtracted, averaged
#(bgsub_norm_averages, bgsub_norm_stds) = averages(bgsub_norm_wells, layout)

# First timepoint shifted to 0 (better for fitting)
#reset_bgsub_means = reset_first_timepoint_to_zero(bgsub_norm_averages)
"""Timecourses normalized, BG-subtracted, averaged, then with first point
shifted to t = 0."""
#reset_bgsub_sds = reset_first_timepoint_to_zero(bgsub_norm_stds)

lipo_conc_conv_factor = 15.502 # 1 mg/ml ~= 15.502 nM liposomes
bg_tc = bgsub_averages['Bax 185 nM, Lipos 0 mg/ml'][VALUE]
bg_time = bgsub_averages['Bax 185 nM, Lipos 0 mg/ml'][TIME]
data_to_fit = []
lipo_concs_to_fit = []
lipo_mgs_to_fit = [1., 0.5, 0.25, 0.125, 0.063, 0.031, 0.016, 0.008, 0]
for conc_name in bgsub_averages.keys():
    lipo_mg = float(conc_name.split()[4])
    if not lipo_mg in lipo_mgs_to_fit:
        continue
    lipo_concs_to_fit.append(lipo_mg * lipo_conc_conv_factor)
    t = bgsub_averages[conc_name][TIME]
    v = bgsub_averages[conc_name][VALUE]
    v_bg = v / bg_tc
    data_to_fit.append(v_bg)

def plot_data():
    """Plots the data and various transformations of it."""
    ion()

    # Timecourse wells
    figure()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    # Lipo bg wells
    figure()
    plot_all(lipo_bg_wells)
    title("Lipo background wells")

    # Background-subtracted wells
    figure()
    plot_all(bgsub_wells)
    title("Background-subtracted wells")

def plot_lipo_background(wells, layout):
    """Takes the lipo background well timecourses and plots the
    fluorescence values of the initial points as a function of lipo
    concentration."""
    init_vals = np.zeros(len(layout.keys()))
    conc_list = np.zeros(len(layout.keys()))
    for i, cond_name in enumerate(layout):
        conc = float(cond_name.split(' ')[1])
        well_name = layout[cond_name][0]
        init_val = wells[well_name][VALUE][0]
        init_vals[i] = init_val
        conc_list[i] = conc

    # Fit the liposome background to a line
    print "Fitting liposome background"
    m = fitting.Parameter(1.)
    b = fitting.Parameter(1.)
    def linear(s):
        return m()*s + b()
    result = fitting.fit(linear, [m, b], init_vals, conc_list)

    plt.figure()
    plt.plot(conc_list, init_vals, marker='o', linestyle='', color='b')
    plt.plot(conc_list, linear(conc_list), color='r')
    plt.xlabel('Liposome concentration (mg/ml)')
    plt.ylabel('RFU')
    plt.title('Liposome background fluorescence at t=0')
    print "m: %f" % m()
    print "b: %f" % b()
    return [m(), b()]

def plot_two_exp_fits():
    data = bgsub_wells
    fmax_arr = np.zeros((1, len(layout.keys())))
    k1_arr = np.zeros((1, len(layout.keys())))
    k2_arr = np.zeros((1, len(layout.keys())))
    conc_list = np.zeros(len(layout.keys()))

    for i, conc_str in enumerate(layout.keys()):
        plt.figure()
        wells = layout[conc_str]
        well_conc = float(conc_str.split(' ')[1])
        conc_list[i] = well_conc

        for j, well in enumerate(wells):
            well_data = data[well]
            time = np.array(well_data[TIME])
            y = np.array(well_data[VALUE])

            fit = TwoExp()
            (k1, fmax, k2) = fit.fit_timecourse(time, y)

            fmax_arr[j, i] = fmax
            k1_arr[j, i] = k1
            k2_arr[j, i] = k2

            plt.plot(time, y, 'b')
            plt.plot(time, fit.fit_func(time, (k1, fmax, k2)), 'r')
            plt.title(well_conc)
            plt.xticks([0, 2000, 4000, 6000, 8000])
            plt.ylabel('% Release')
            plt.xlabel('Time (sec)')

    plt.show()
    return (fmax_arr, k1_arr, k2_arr, conc_list)

def plot_k1_curve(k1_arr, conc_list):
    k1_means = np.mean(k1_arr, axis=0)
    k1_sds = np.std(k1_arr, axis=0)

    plt.figure()
    plt.errorbar(conc_list, k1_means, yerr= k1_sds / np.sqrt(3), color='r',
                 linestyle='', linewidth=2)
    plt.title('$k_1$')
    plt.xlabel('[Bax] (nM)')
    plt.ylabel('$k_1\ (\mathrm{sec}^{-1})$')
    """
    # Fit with exponential-linear curve
    vi = fitting.Parameter(0.05)
    vf = fitting.Parameter(0.025)
    tau = fitting.Parameter(0.02)
    # Define fitting function
    def biphasic(t):
        return (vf()*t) + ( (vi() - vf()) *
                            ((1 - np.exp(-tau()*t))/tau()) )
    fitting.fit(biphasic, [vi, vf, tau], k1_means, conc_list)
    plt.plot(conc_list, biphasic(conc_list), 'r', linewidth=2)

    # Fit with "double-binding curve"
    ksat = fitting.Parameter(0.0005)
    knonsat = fitting.Parameter(20)
    R0 = fitting.Parameter(20)
    def double_binding(t):
        return (R0() * t)/(ksat() + t) + (knonsat()*t)
    fitting.fit(double_binding, [R0, ksat, knonsat], k1_means, conc_list)
    plt.plot(conc_list, double_binding(conc_list), 'g', linewidth=2)

    plt.text(400, 0.00025,
            r'$\frac{R_0 Bax_0}{K_{sat} + Bax_0} + K_{nonsat} Bax_0$')
    plt.text(400, 0.00020, r'$R_0$ = %.3g nM' % R0())
    plt.text(400, 0.00015, r'$K_{sat}$ = %.3g nM' % ksat())
    plt.text(400, 0.00010, r'$K_{nonsat}$ = %.3g nM' % knonsat())

    plt.text(35, 0.00054,
            r'$v_f + (v_i - v_f) \left(\frac{1 - e^{-\tau t}}{\tau}\right)$')
    plt.text(35, 0.00049, r'$v_i$ = %.3g sec$^{-1}$ nM$^{-1}$' % vi())
    plt.text(35, 0.00044, r'$v_f$ = %.3g sec$^{-1}$ nM$^{-1}$' % vf())
    plt.text(35, 0.00039, r'$\frac{1}{\tau}$ = %.2f nM' % (1/tau()))
    plt.title('Biphasic fit of $k_1$')
    plt.show()
    """

def plot_k2_curve(k2_arr, conc_list):
    k2_means = np.mean(k2_arr, axis=0)
    k2_sds = np.std(k2_arr, axis=0)
    # Plot k2 on a linear scale
    plt.figure()
    plt.errorbar(conc_list, k2_means, yerr=k2_sds / np.sqrt(3))
    plt.title('$k_2$')

def plot_fmax_curve(fmax_arr, conc_list):
    fmax_means = np.mean(fmax_arr, axis=0)
    fmax_sds = np.std(fmax_arr, axis=0)

    # Plot of Fmax
    plt.figure()
    plt.ylabel(r'$F_{max}$ value (% Release)')
    plt.xlabel('[Bax] (nM)')
    plt.title(r'$F_{max}$')

    # Try fitting the fmax values to a hill function
    hill_vmax = fitting.Parameter(1)
    kd = fitting.Parameter(100)
    def hill(s):
        return ((hill_vmax() * s) / (kd() + s))
    fitting.fit(hill, [kd], fmax_means[1:], conc_list[1:])
    plt.plot(conc_list, hill(conc_list), 'r', linewidth=2, label='Hill fit')

    plt.text(35, 0.93, '$K_D$ = %.2f' % kd(), fontsize=16)
    #plt.text(35, 0.99, '$V_{max}$ = %.2f' % hill_vmax(), fontsize=16)

    # Try fitting the fmax values to an exponential function
    exp_fmax = fitting.Parameter(1.0)
    exp_k = fitting.Parameter(0.01)
    def exp_func(s):
        return exp_fmax() * (1 - np.exp(-exp_k()*s))
    fitting.fit(exp_func, [exp_fmax, exp_k], fmax_means[1:], conc_list[1:])
    plt.plot(conc_list, exp_func(conc_list), 'g', linewidth=2,
             label='Exp fit')

    # Plot the data
    plt.errorbar(conc_list[1:], fmax_means[1:],
                 yerr=fmax_sds[1:] / np.sqrt(3),
                 label='Data', linestyle='', linewidth=2)

    plt.legend(loc='lower right')
    plt.show()

def mpi_fit():
    bg_time = bgsub_averages['Bax 185 nM, Lipos 0 mg/ml'][TIME]
    bg_tc = bgsub_averages['Bax 185 nM, Lipos 0 mg/ml'][VALUE]
    lipo_concs = []
    data = []
    lipo_concs_to_fit = [1., 0.5, 0.25, 0.125, 0.063, 0.031]
    for conc_name in bgsub_averages.keys():
        lipo_conc = float(conc_name.split()[4])
        if not lipo_conc in lipo_concs_to_fit:
            continue
        lipo_concs.append(lipo_conc * 15.502)
        t = bgsub_averages[conc_name][TIME]
        v = bgsub_averages[conc_name][VALUE]
        v_bg = v / bg_tc
        data.append(v_bg)

    (gf, sampler) = fit_with_2conf_mc(bg_time, data, lipo_concs)
    #gf = fit_with_2conf(bg_time, data, lipo_concs)
    #gf = fit_with_2conf_rev(bg_time, data, lipo_concs)
    #gf = fit_with_3conf(bg_time, data, lipo_concs)
    #gf = fit_with_simple_3conf(bg_time, data, lipo_concs)
    return (gf, sampler)

def fit_with_simple_3conf(time, data, lipo_concs):
    params_dict = {
                   'c1_scaling': 2.,
                   'c2_scaling': 4.85,
                   'c0_to_c1_k': 1e-3,
                   'c1_to_c2_k': 1e-3,
                  }
    bd = multiconf.Builder(params_dict=params_dict)
    bd.build_model_multiconf(3, 1, normalized_data=True, reversible=False)
    bd.global_params = (bd['c1_scaling'],
                        bd['c2_scaling'],
                        bd['c1_to_c2_k'])
    bd.local_params = [
                        bd['c0_to_c1_k'],
            ]
    params = {}
    gf = fitting.GlobalFit(bd, time, data, params, 'NBD')
    gf.fit()
    return gf

def fit_with_3conf(time, data, lipo_concs):
    params_dict = {
                   'c1_scaling': 1.,
                   'c2_scaling': 4.8,
                   'Bax_transloc_kf': 1e-4,
                   'Bax_transloc_kr': 1e-1,
                   'basal_Bax_kf': 9e-1,
                   'pore_formation_rate_k': 3e-3,
                  }
    bd = one_cpt.Builder(params_dict=params_dict)
    bd.build_model_nbd_3_conf()
    bd.global_params = (bd['c1_scaling'],
                        bd['c2_scaling'],
                        bd['Bax_transloc_kf'],
                        bd['Bax_transloc_kr'],
                        #bd['basal_Bax_kf'],
                        bd['pore_formation_rate_k'])
    bd.local_params = [bd['basal_Bax_kf']]
    params = {'Vesicles_0': lipo_concs}
    gf = fitting.GlobalFit(bd, time, data, params, 'NBD')
    gf.fit()
    return gf

def fit_with_2conf(time, data, lipo_concs):
    params_dict = {'c1_scaling': 4.85,
                   'basal_Bax_kf': 1.73e-3,
                   'Bax_transloc_kr': 1.4835e-1,
                   'Bax_transloc_kf': 1e-2}
    bd = one_cpt.Builder(params_dict=params_dict)
    bd.build_model_nbd_2_conf()
    bd.global_params = (bd['c1_scaling'],
                        bd['Bax_transloc_kf'],
                        bd['Bax_transloc_kr'],
                        bd['basal_Bax_kf'])
    bd.local_params = []
    params = {'Vesicles_0': lipo_concs}
    gf = fitting.GlobalFit(bd, time, data, params, 'NBD')
    gf.fit()
    return gf

def fit_with_2conf_rev(time, data, lipo_concs):
    params_dict = {'c1_scaling': 4.85,
                   'basal_Bax_kf': 1.73e-3,
                   'Bax_transloc_kr': 1.4835e-1,
                   'Bax_transloc_kf': 1e-2,
                   'iBax_reverse_k': 2e-4}
    bd = one_cpt.Builder(params_dict=params_dict)
    bd.build_model_nbd_2_conf_rev()
    bd.global_params = (bd['c1_scaling'],
                        bd['Bax_transloc_kf'],
                        bd['Bax_transloc_kr'],
                        bd['basal_Bax_kf'],
                        bd['iBax_reverse_k'])
    bd.local_params = []
    params = {'Vesicles_0': lipo_concs}
    gf = fitting.GlobalFit(bd, time, data, params, 'NBD')
    gf.fit()
    return gf

def plot_emcee_fits(gf, sampler):
    set_fig_params_for_publication()
    fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
    plt.ylabel('$F/F_0$')
    plt.xlabel(r'Time (sec $\times 10^3$)')
    plt.ylim([0.7, 5.2])
    plt.xlim([0, bg_time[-1] + 500])
    ax = plt.gca()
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.yaxis.set_tick_params(direction='out')
    ax.xaxis.set_tick_params(direction='out')
    ax.set_xticks(np.linspace(0, 1e4, 6))
    ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
    plt.subplots_adjust(bottom=0.24, left=0.21)
    for data in data_to_fit:
        plt.plot(bg_time, data, 'k', linewidth=1)
    # Plot the final point
    gf.plot_func(sampler.flatchain[0,-1,:])

if __name__ == '__main__':
    import triangle

    if len(sys.argv) == 3:
        mcmc_path = sys.argv[1]
        output_base = sys.argv[2]
        # Load the sampler data
        (gf, sampler) = pickle.load(open(mcmc_path))
        # Plot the best fit vs. the data
        plot_emcee_fits(gf, sampler)
        plt.savefig('%s_fits.pdf' % output_base)
        # Triangle plot for parameters
        plt.figure()
        triangle.corner(sampler.flatchain[0])
        plt.savefig('%s_tri.pdf' % output_base)
        sys.exit()

    (gf2, sampler2) = pickle.load(open('../pt/pt_140318_nbd_2_conf_4.pck'))
    plot_emcee_fits(gf2, sampler2)
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(sampler2._lnprob[0,:,:].T, alpha=0.1)
    plt.title('0th chain')
    plt.subplot(2, 2, 2)
    plt.plot(sampler2._lnprob[1,:,:].T, alpha=0.1)
    plt.title('1th chain')
    plt.subplot(2, 2, 3)
    plt.plot(sampler2._lnprob[2,:,:].T, alpha=0.1)
    plt.title('2th chain')
    plt.subplot(2, 2, 4)
    plt.plot(sampler2._lnprob[3,:,:].T, alpha=0.1)
    plt.title('3th chain')

    (gf_2r, sampler_2r) = pickle.load(open('../pt/pt_140318_nbd_2_conf_rev_3.pck'))
    plot_emcee_fits(gf_2r, sampler_2r)
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(sampler_2r._lnprob[0,:,:].T, alpha=0.1)
    plt.title('0th chain')
    plt.subplot(2, 2, 2)
    plt.plot(sampler_2r._lnprob[1,:,:].T, alpha=0.1)
    plt.title('1th chain')
    plt.subplot(2, 2, 3)
    plt.plot(sampler_2r._lnprob[2,:,:].T, alpha=0.1)
    plt.title('2th chain')
    plt.subplot(2, 2, 4)
    plt.plot(sampler_2r._lnprob[3,:,:].T, alpha=0.1)
    plt.title('3th chain')

    (gf3, sampler3) = pickle.load(open('../pt/pt_140318_nbd_3_conf_2.pck'))
    plot_emcee_fits(gf3, sampler3)
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(sampler3._lnprob[0,:,:].T, alpha=0.1)
    plt.title('0th chain')
    plt.subplot(2, 2, 2)
    plt.plot(sampler3._lnprob[1,:,:].T, alpha=0.1)
    plt.title('1th chain')
    plt.subplot(2, 2, 3)
    plt.plot(sampler3._lnprob[2,:,:].T, alpha=0.1)
    plt.title('2th chain')
    plt.subplot(2, 2, 4)
    plt.plot(sampler3._lnprob[3,:,:].T, alpha=0.1)
    plt.title('3th chain')

    sys.exit()

    import triangle
    import pickle

    #(gf, sampler) = plot_timecourse_figure()
    (gf, sampler) = mpi_fit()
    chain = sampler.flatchain
    #fig = triangle.corner(chain)
    #fig.savefig("triangle0.png")

    with open('140318fit_pt.pck', 'w') as f:
        pickle.dump((gf, chain), f)

    #plt.figure()
    #plt.hist(chain[:,0], bins=20)
    sys.exit()

    #plot_data()
    #plot_lipo_background(lipo_bg_wells, lipo_bg_layout)

    # Try fitting the high conc trajectory
    y = bgsub_wells['A3'][VALUE]
    t = bgsub_wells['A3'][TIME]
    b = fitting.Parameter(np.mean(y[0:2]))
    fmax = fitting.Parameter(25.)
    k1 = fitting.Parameter(np.log(2)/2000.)
    k2 = fitting.Parameter(1e-3)
    k_bleach = fitting.Parameter(1.17e-5)
    def exp_func(t):
        return (b() + fmax()*(1 - np.exp(-k1()*(1 - np.exp(-k2()*t))*t))) * \
                np.exp(-k_bleach()*t)

    # One-parameter exp
    #def exp_func(t):
    #    return (b() + fmax()*(1 - np.exp(-k1()*t))) * \
    #            np.exp(-k_bleach()*t)

    fitting.fit(exp_func, [fmax, k1, k2], y, t)

    plt.figure()
    plt.plot(t, y, marker='o', linestyle='')
    plt.plot(t, exp_func(t), color='r', linewidth='2')

    print "Fmax: %f" % ((fmax() + b()) / b())

    sys.exit()
    #(fmax_arr, k1_arr, k2_arr, conc_list) = plot_two_exp_fits()
    #plot_fmax_curve(fmax_arr, conc_list)

    from tbidbaxlipo.plots import titration_fits
    """
    fit = TwoExp()
    pore_df = to_dataframe(pores, bgsub_norm_stds)
    fit.plot_fits_from_dataframe(pore_df)
    p = fit.fit_from_dataframe(pore_df)


    # Multiply Fmax values by molar concentration of liposomes
    concs = np.array(pore_df.columns.values, dtype='float')[1:-1]
    concentration_of_pores = p[1][1:-1] * 0.775
    plt.figure()
    plt.plot(concs, concentration_of_pores)

    # Fit a straight line to the concentrations
    from tbidbaxlipo.util import fitting
    m = fitting.Parameter(0.025)
    b = fitting.Parameter(0)
    def linear(x):
        return m()*x + b()
    fitting.fit(linear, [m, b], concentration_of_pores, concs)
    plt.plot(concs, linear(concs))
    """

    df = to_dataframe(norm_averages, norm_stds)
    concs = np.array(df.columns.values, dtype='float')
    fit = titration_fits.TwoExp()

    # Get background average
    background_time = norm_averages['Bax 0 nM'][TIME]
    background = norm_averages['Bax 0 nM'][VALUE]

    bg_rate = fit.fit_timecourse(background_time, background)

    plt.ion()
    plt.plot(background_time, background)
    plt.plot(background_time, fit.fit_func(background_time, bg_rate))
    t = np.linspace(0, 100000, 1000)
    plt.plot(t, fit.fit_func(t, bg_rate))
    import pdb; pdb.set_trace()

    fit = titration_fits.TwoExp()
    #fit.plot_fits_from_dataframe(subset_df)
    #p = fit.fit_from_dataframe(subset_df)
    fit.plot_fits_from_dataframe(df)
    p = fit.fit_from_dataframe(df)

    # With fitting of bg
    #print "bg_rate %f" % bg_rate
    fit = titration_fits.TwoExpWithBackground(bg_rate)
    #fit.plot_fits_from_dataframe(subset_df)
    #p = fit.fit_from_dataframe(subset_df)
    fit.plot_fits_from_dataframe(df)
    p = fit.fit_from_dataframe(df)
