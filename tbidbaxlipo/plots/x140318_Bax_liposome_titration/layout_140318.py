import sys
import pickle
import numpy as np
from matplotlib import pyplot as plt
from preprocess_data import timecourse_wells, lipo_bg_wells, bgsub_wells, \
                            layout, lipo_bg_layout, bax_lipo_layout, \
                            data_to_fit, bg_time, lipo_concs_to_fit
from tbidbaxlipo.util.plate_assay import plot_all, TIME, VALUE
from tbidbaxlipo.plots.titration_fits import TwoExp, OneExpFmax
from tbidbaxlipo.util import fitting, set_fig_params_for_publication, \
                             emcee_fit, format_axis
from tbidbaxlipo.models.nbd import multiconf
from tbidbaxlipo.models import one_cpt

def plot_data():
    """Plots the data and various transformations of it."""
    # Timecourse wells
    plt.figure()
    plot_all(timecourse_wells)
    plt.title("Raw timecourses")

    # Lipo bg wells
    plt.figure()
    plot_all(lipo_bg_wells)
    plt.title("Lipo background wells")

    # Background-subtracted wells
    plt.figure()
    plot_all(bgsub_wells)
    plt.title("Background-subtracted wells")

def plot_lipo_background(wells, layout):
    """Takes the lipo background well timecourses and plots the
    fluorescence values of the initial points as a function of lipo
    concentration. Shows that the liposome background fluorescence
    is strictly linear.
    """
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
    ax = plt.gca()
    ax.set_xscale('log')
    print "m: %f" % m()
    print "b: %f" % b()
    return [m(), b()]

def plot_exp_fits(time, data, concs, plot_fits=True):
    fmax_arr = np.zeros(len(concs))
    k_arr = np.zeros(len(concs))
    for i, conc in enumerate(concs):
        y = data[i]
        fmax = fitting.Parameter(3.)
        k = fitting.Parameter(5e-4)
        def fit_func(t):
            return 1 + fmax() * (1 - np.exp(-k()*t))
        fitting.fit(fit_func, [k, fmax], y, time)

        fmax_arr[i] = fmax()
        k_arr[i] = k()

        if plot_fits:
            plt.figure()
            plt.plot(time, y, 'b')
            plt.plot(time, fit_func(time), 'r')
            plt.title('%f nM liposomes' % conc)
            plt.xticks([0, 2000, 4000, 6000, 8000])
            plt.ylabel('$F/F_0$')
            plt.xlabel('Time (sec)')
            plt.show()

    return (fmax_arr, k_arr)

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
    ax.set_xticks(np.linspace(0, 1e4, 6))
    ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
    plt.subplots_adjust(bottom=0.24, left=0.21)
    for data in data_to_fit:
        plt.plot(bg_time, data, 'k', linewidth=1)
    # Plot the final point
    gf.plot_func(sampler.flatchain[0,-1,:])
    format_axis(ax)

if __name__ == '__main__':
    import triangle
    plt.ion()

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

