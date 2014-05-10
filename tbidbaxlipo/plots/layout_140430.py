import sys
import os
import tbidbaxlipo.data
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.models.nbd.multiconf import Builder
from tbidbaxlipo.util import fitting

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '140430_cBid85CDAC_BaxNBD_FRET.csv'))

data_arr = np.loadtxt(timecourse_file, delimiter=',')


data = {
        'c36_nbd_time': data_arr[:,0],
        'c36_nbd': data_arr[:,1],
        'c36_fret_time': data_arr[:,2],
        'c36_fret': data_arr[:,3],
        'c68_nbd_time': data_arr[:,4],
        'c68_nbd': data_arr[:,5],
        'c68_fret_time': data_arr[:,6],
        'c68_fret': data_arr[:,7],
        'c126_nbd_time': data_arr[:,8],
        'c126_nbd': data_arr[:,9],
        'c126_fret_time': data_arr[:,10],
        'c126_fret': data_arr[:,11]}

def plot_fret():
    plt.figure('FRET')
    plt.plot(data['c36_fret_time'], data['c36_fret'], label='c36')
    plt.plot(data['c68_fret_time'], data['c68_fret'], label='c68')
    plt.plot(data['c126_fret_time'], data['c126_fret'], label='c126')
    plt.xlabel('Time (sec)')
    plt.ylabel('% FRET')
    plt.title('DAC-Bid/NBD-Bax FRET')
    plt.legend(loc='right')

def plot_nbd():
    plt.figure('NBD')
    plt.plot(data['c36_nbd_time'], data['c36_nbd'], label='c36')
    plt.plot(data['c68_nbd_time'], data['c68_nbd'], label='c68')
    plt.plot(data['c126_nbd_time'], data['c126_nbd'], label='c126')
    plt.xlabel('Time (sec)')
    plt.ylabel('$F/F_0$')
    plt.title('NBD Bax Fluorescence')
    plt.legend(loc='right')

def plot_fits():
    params_dict = {'c0_to_c1_k': 2e-3,
                   'c1_scaling': 0.4,
                   'c1_to_c2_k': 1e-3,
                   'c2_scaling': 0.6,
                   'c1_to_c2_k': 1e-3,
                   'c3_scaling': 0.5}

    builder = Builder(params_dict=params_dict)
    builder.build_model_multiconf(5, data['c68_fret'][0],
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

    pysb_fit = fitting.fit_pysb_builder(builder, 'NBD', data['c68_fret_time'],
                                        data['c68_fret'])

    plt.figure()
    plt.plot(data['c68_fret_time'], data['c68_fret'], linestyle='', marker='.')
    plt.plot(data['c68_fret_time'], pysb_fit.ypred)
    plt.xlabel('Time (sec)')

if __name__ == '__main__':
    plt.ion()
    plot_fret()
    plot_nbd()
    plot_fits()

