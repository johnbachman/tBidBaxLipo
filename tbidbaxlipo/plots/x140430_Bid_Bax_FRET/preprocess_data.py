import os
import sys
import numpy as np
import tbidbaxlipo.data
from tbidbaxlipo.util.calculate_error_variance import calc_err_var

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

def nbd_fret_err(nbd_data, fret_data, plot=False, plot_title=None):
    """Utility function for calculating NBD and FRET data error."""
    # Plot titles
    if plot_title:
        nbd_plot_title = plot_title + ' NBD'
        fret_plot_title = plot_title + ' Bid/Bax FRET'
    else:
        nbd_plot_title = None
        fret_plot_title = None
    # Calculate NBD err
    nbd_residuals, _ = calc_err_var(nbd_data, last_n_pts=80, plot=plot,
                                    plot_title=nbd_plot_title)
    nbd_err = np.std(nbd_residuals, ddof=1)
    # Calculate FRET err
    fret_residuals, _ = calc_err_var(fret_data, last_n_pts=80, plot=plot,
                                     plot_title=fret_plot_title)
    fret_err = np.std(fret_residuals, ddof=1)
    return (nbd_err, fret_err)

time_36 = data_dict['c36_fret_time'] # Same as NBD time
data_36 = np.array([[data_dict['c36_nbd'], data_dict['c36_fret']]])
c36_err = nbd_fret_err(data_dict['c36_nbd'], data_dict['c36_fret'])
data_sigma_36 = np.array([c36_err])

time_68 = data_dict['c68_fret_time'] # Same as NBD time
data_68 = np.array([[data_dict['c68_nbd'], data_dict['c68_fret']]])
c68_err = nbd_fret_err(data_dict['c68_nbd'], data_dict['c68_fret'])
data_sigma_68 = np.array([c68_err])

time_126 = data_dict['c126_fret_time'] # Same as NBD time
data_126 = np.array([[data_dict['c126_nbd'], data_dict['c126_fret']]])
c126_err = nbd_fret_err(data_dict['c126_nbd'], data_dict['c126_fret'])
data_sigma_126 = np.array([c126_err])


