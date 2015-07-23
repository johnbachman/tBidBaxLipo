from tbidbaxlipo.pt import run_pt

fit_dict = {
    'data': {
        'module': 'tbidbaxlipo.plots.x140318_Bax_liposome_titration.preprocess_data',
        'data_var': 'data_to_fit',
        'data_sigma_var': 'data_sigma',
        'initial_condition_var': 'lipo_concs_to_fit',
        'time_var': 'bg_time',
    },
    'model': {
        'baxtranslocation': 1,
        'activation': 1,
        'nbd': 1,
    },
    'model_observable': 'NBD',
    'global_initial_conditions': {
        'tBid_0': 0,
        'c0_scaling': 1.0,
        'Bax_0': 185.0,
    },
    'local_initial_condition': 'Vesicles_0',
    'global_params': 'all',
    'local_params': [],
    'ntemps': 3,
    'highest_temp': -2,
    'nwalkers': 50,
    'nburnin': 1,
    'nsample': 1,
    'thin': 1,
}

def test_run_pt():
    run_pt.run(fit_dict)
