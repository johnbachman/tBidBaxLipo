import os
import sys
import numpy as np
import tbidbaxlipo.data

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '140430_cBid85CDAC_BaxNBD_FRET.csv'))

data_arr = np.loadtxt(timecourse_file, delimiter=',')

data_dict = {
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

# First axis of data array: conditions (e.g. concentrations)
# Second axis: observables
# Third axis: time

time_36 = data_dict['c36_fret_time'] # Same as NBD time
data_36 = np.array([[data_dict['c36_nbd'], data_dict['c36_fret']]])

time_68 = data_dict['c68_fret_time'] # Same as NBD time
data_68 = np.array([[data_dict['c68_nbd'], data_dict['c68_fret']]])

time_126 = data_dict['c126_fret_time'] # Same as NBD time
data_126 = np.array([[data_dict['c126_nbd'], data_dict['c126_fret']]])


