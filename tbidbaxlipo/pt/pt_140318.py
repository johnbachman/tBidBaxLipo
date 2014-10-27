import sys
import pickle
from tbidbaxlipo.models import one_cpt
from tbidbaxlipo.util import emcee_fit
from tbidbaxlipo.plots.layout_140318 import data_to_fit, lipo_concs_to_fit, \
                                            bg_time

def fit(model_name, time, data, lipo_concs):
    params_dict = {'Bax_0': 185.,
                   'c1_scaling': 4.85,
                   'basal_Bax_kf': 1.73e-3,
                   'Bax_transloc_kr': 1.4835e-1,
                   'Bax_transloc_kf': 1e-2}
    bd = one_cpt.Builder(params_dict=params_dict)
    build_model_cmd = 'bd.build_model_%s()' % model_name
    eval(build_model_cmd)
    bd.global_params = (bd['c1_scaling'],
                        bd['Bax_transloc_kf'],
                        bd['Bax_transloc_kr'],
                        bd['basal_Bax_kf'])
    bd.local_params = []
    params = {'Vesicles_0': lipo_concs}
    gf = emcee_fit.GlobalFit(bd, time, data, params, 'NBD')
    sampler = emcee_fit.pt_mpi_sample(gf, 25, 100, 5, 10)
    return (gf, sampler)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Please include the model name as an argument."
        sys.exit()

    model_name = sys.argv[1]
    (gf, sampler) = fit(model_name, bg_time, data_to_fit, lipo_concs_to_fit)

    with open('pt_140318_%s.pck' % model_name, 'w') as f:
        pickle.dump((gf, sampler), f)


