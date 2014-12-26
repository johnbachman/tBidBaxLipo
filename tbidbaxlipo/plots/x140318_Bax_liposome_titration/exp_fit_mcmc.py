import numpy as np
from matplotlib import pyplot as plt

from preprocess_data import data_to_fit, bg_time, lipo_concs_to_fit
from tbidbaxlipo.util import set_fig_params_for_publication, \
                             emcee_fit, format_axis
from pysb import *
from tbidbaxlipo.models.nbd.multiconf import Builder

bd = Builder()
bd.build_model_multiconf(num_confs=2, c0_scaling=1, normalized_data=True,
                        reversible=False)

def fit_with_2conf(time, data, lipo_concs):
    y = [data[-2]]
    bd.global_params = [bd['c1_scaling'], bd['c0_to_c1_k']]
    bd.local_params = []
    params = {}
    gf = emcee_fit.GlobalFit(bd, time, y, params, 'NBD')
    sampler = emcee_fit.ens_sample(gf, 100, 20, 20)
    return (gf, sampler)

if __name__ == '__main__':
    (gf, sampler) = fit_with_2conf(bg_time, data_to_fit, lipo_concs_to_fit)

