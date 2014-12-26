import sys
import pickle
import numpy as np
from matplotlib import pyplot as plt
from preprocess_data import timecourse_wells, lipo_bg_wells, bgsub_wells, \
                            layout, lipo_bg_layout, bax_lipo_layout, \
                            data_to_fit, bg_time, lipo_concs_to_fit
from tbidbaxlipo.util.plate_assay import plot_all, TIME, VALUE
from tbidbaxlipo.util import fitting, set_fig_params_for_publication, \
                             emcee_fit, format_axis
from tbidbaxlipo.models.nbd import multiconf
from tbidbaxlipo.models import one_cpt

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
