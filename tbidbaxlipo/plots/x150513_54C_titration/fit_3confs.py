from tbidbaxlipo.plots.x150513_54C_titration.preprocess_data import data_dict
from tbidbaxlipo.util import fitting
from tbidbaxlipo.models.nbd.multiconf import Builder
from matplotlib import pyplot as plt
import numpy as np

# -- Model fits --
params_dict = {'c1_to_c2_k': 1e-4, 'c1_scaling': 2,
               'c0_to_c1_k': 2e-3}

def fit_3confs(t, y, do_plot=True):
    # Set up the NBD/FRET model
    builder = Builder(params_dict=params_dict)
    builder.build_model_multiconf(3, y[0], normalized_data=True)
    # Add some initial guesses
    # Guess that the scaling value for the final conformation is close
    # to the final data value
    builder.model.parameters['c2_scaling'].value = y[-1]
    # Guess that the scaling value for the intermediate conformation
    # is close to the value at ~200 sec
    c1_timescale_seconds = 300
    c1_timescale_index = np.where(t > c1_timescale_seconds)[0].min()
    builder.model.parameters['c1_scaling'].value = \
                                            y[c1_timescale_index]
    # Rough guesses for the timescales of the first and second
    # transitions
    builder.model.parameters['c0_to_c1_k'].value = 0.025
    builder.model.parameters['c1_to_c2_k'].value = 5e-3
    print builder.model.parameters
    # Do the fit
    pysb_fit = fitting.fit_pysb_builder(builder, 'NBD', t, y,
                                        log_transform=True)
    import ipdb; ipdb.set_trace()
    if do_plot:
        plt.figure('fits')
        plt.plot(t, y, linestyle='', marker='.', markersize=10)
        plt.plot(t, pysb_fit.ypred, color='k')
        plt.xlabel('Time (sec)')
        plt.ylabel('NBD F/F0')
    return pysb_fit

def local_fits():
    fit62 = fit_3confs(data_dict['nt_62'], data_dict['ny_62'])
    fit125 = fit_3confs(data_dict['nt_125'], data_dict['ny_125'])
    fit250 = fit_3confs(data_dict['nt_250'], data_dict['ny_250'])
    fit500 = fit_3confs(data_dict['nt_500'], data_dict['ny_500'])
    return [fit62, fit125, fit250, fit500]

def global_fluorescence_det_fits(data_dict):
    t = data_dict['nt_62']
    data_list = [data_dict[key] for key in
                 ['ny_62', 'ny_125', 'ny_250', 'ny_500']]
    # Set up the NBD/FRET model
    builder = Builder(params_dict=params_dict)
    #builder.build_model_multiconf(3, 1., normalized_data=True, reversible=True)
    builder.build_model_3confs_c2_dimer(1., normalized_data=True, reversible=True)
    # Guess a c1_scaling value drawn from the 62.5 nM data
    c1_timescale_seconds = 300
    c1_timescale_index = np.where(t > c1_timescale_seconds)[0].min()
    builder.model.parameters['c1_scaling'].value = \
                            data_dict['ny_62'][c1_timescale_index]
    # Guess a c2_scaling value drawn from the 500 nM data
    builder.model.parameters['c2_scaling'].value = data_dict['ny_500'][-1]
    # Rough guesses for the timescales of the first and second
    # transitions
    builder.model.parameters['c0_to_c1_k'].value = 0.025
    builder.model.parameters['c1_to_c2_k'].value = 5e-3
    # Set the list of global and local params
    global_param_names = ['c0_scaling', 'c1_scaling', 'c2_scaling']
    local_param_names = ['c0_to_c1_k', 'c1_to_c2_k', 'c1_to_c0_k',
                         'c2_to_c1_k']
    builder.global_params = [builder.model.parameters[name]
                             for name in global_param_names]
    builder.local_params = [builder.model.parameters[name]
                            for name in local_param_names]

    # The list of initial conditions
    params = {'Bax_0': [62.5, 125, 250, 500]}
    # Build the global fit object
    gf = fitting.GlobalFit(builder, t, data_list, params, 'NBD',
                           obs_type='Expression')
    gf.fit()
    return gf

if __name__ == '__main__':
    plt.ion()
    gf = global_fluorescence_fits(data_dict)
    plt.figure()
    t = data_dict['nt_62']
    data_list = [data_dict[key] for key in
                 ['ny_62', 'ny_125', 'ny_250', 'ny_500']]
    for y in data_list:
        plt.plot(t, y, linestyle='', marker='.', markersize=10, color='b')
    gf.plot_func()


