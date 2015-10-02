import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.models.nbd.multiconf import Builder
from tbidbaxlipo.util import fitting
from preprocess_data import data_dict, nbd_fret_err

def plot_fret():
    plt.figure('FRET')
    plt.plot(data_dict['c36_fret_time'], data_dict['c36_fret'], label='c36')
    plt.plot(data_dict['c68_fret_time'], data_dict['c68_fret'], label='c68')
    plt.plot(data_dict['c126_fret_time'], data_dict['c126_fret'], label='c126')
    plt.xlabel('Time (sec)')
    plt.ylabel('% FRET')
    plt.title('DAC-Bid/NBD-Bax FRET')
    plt.legend(loc='right')

def plot_nbd():
    plt.figure('NBD')
    plt.plot(data_dict['c36_nbd_time'], data_dict['c36_nbd'], label='c36')
    plt.plot(data_dict['c68_nbd_time'], data_dict['c68_nbd'], label='c68')
    plt.plot(data_dict['c126_nbd_time'], data_dict['c126_nbd'], label='c126')
    plt.xlabel('Time (sec)')
    plt.ylabel('$F/F_0$')
    plt.title('NBD Bax Fluorescence')
    plt.legend(loc='right')

def plot_error_estimates():
    nbd_fret_err(data_dict['c36_nbd'], data_dict['c36_fret'], plot=True,
                 plot_title='Est. error for NBD-C36-Bax, ')
    nbd_fret_err(data_dict['c68_nbd'], data_dict['c68_fret'], plot=True,
                 plot_title='Est. error for NBD-C68-Bax, ')
    nbd_fret_err(data_dict['c126_nbd'], data_dict['c126_fret'], plot=True,
                 plot_title='Est. error for NBD-C126-Bax, ')

def plot_fits():
    params_dict = {'c0_to_c1_k': 2e-2,
                   'c1_scaling': 20.,
                   'c1_to_c2_k': 1e-4,
                   'c2_scaling': 6.,
                   'c1_to_c2_k': 1e-3,
                   'c3_scaling': 0.5}

    builder = Builder(params_dict=params_dict)
    builder.build_model_multiconf(3, data_dict['c126_nbd'][0],
                                  normalized_data=True,
                                  reversible=False)

    """
    k1 = builder.model.parameters['c0_to_c1_k']
    k2 = builder.model.parameters['c1_to_c2_k']
    k1_index = builder.model.parameters.index(k1)
    k2_index = builder.model.parameters.index(k2)
    k1_est_index = builder.estimate_params.index(k1)
    k2_est_index = builder.estimate_params.index(k2)
    """

    pysb_fit = fitting.fit_pysb_builder(builder, 'NBD', data_dict['c126_nbd_time'],
                                        data_dict['c126_nbd'])

    plt.figure()
    plt.plot(data_dict['c126_nbd_time'], data_dict['c126_nbd'], linestyle='', marker='.')
    plt.plot(data_dict['c126_nbd_time'], pysb_fit.ypred)
    plt.xlabel('Time (sec)')

if __name__ == '__main__':
    plt.ion()
    plot_fits()
    #plot_fret()
    #plot_nbd()
    #plot_error_estimates()
