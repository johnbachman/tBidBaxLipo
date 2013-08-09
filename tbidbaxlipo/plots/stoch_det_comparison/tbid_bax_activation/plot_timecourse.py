from tbidbaxlipo.models import one_cpt, n_cpt, site_cpt, simulation
from pylab import *
import collections
import glob
from scipy.stats import poisson
import pkgutil
import pickle
import numpy as np
import h5py

ion()

def plot_timecourse_comparison(job, data, dt):
    b_one = job.one_cpt_builder()
    [time, one_cpt_obs] = job.run_one_cpt()
    one_cpt_Bax_0 = b_one.model.parameters['Bax_0'].value
    one_cpt_tBid_0 = b_one.model.parameters['tBid_0'].value

    b_n = job.n_cpt_builder()
    n_cpt_means = np.mean(data, axis=0)
    n_cpt_sd = np.std(data, axis=0)
    n_cpt_Bax_0 = b_n.model.parameters['Bax_0'].value
    n_cpt_tBid_0 = b_n.model.parameters['tBid_0'].value

    # Compare deterministic to n_cpt
    figure()
    mBax_se = (n_cpt_sd[dt['mBax']] / n_cpt_Bax_0) / np.sqrt(job.num_sims)
    cBax_se = (n_cpt_sd[dt['cBax']] / n_cpt_Bax_0) / np.sqrt(job.num_sims)
    iBax_se = (n_cpt_sd[dt['iBax']] / n_cpt_Bax_0) / np.sqrt(job.num_sims)
    mtBid_se = (n_cpt_sd[dt['mtBid']] / n_cpt_tBid_0) / np.sqrt(job.num_sims)
    ctBid_se = (n_cpt_sd[dt['ctBid']] / n_cpt_tBid_0) / np.sqrt(job.num_sims)

    n_cpt_time = n_cpt_means[dt['time']]
    errorbar(n_cpt_time,  n_cpt_means[dt['mBax']] / n_cpt_Bax_0, \
             yerr=mBax_se, color='r', label='mBax');
    errorbar(n_cpt_time,  n_cpt_means[dt['cBax']] / n_cpt_Bax_0, \
             yerr=cBax_se, color='b', label='cBax');
    errorbar(n_cpt_time,  n_cpt_means[dt['iBax']] / n_cpt_Bax_0, \
             yerr=iBax_se, color='g', label='iBax');
    errorbar(n_cpt_time,  n_cpt_means[dt['mtBid']] / n_cpt_tBid_0, \
             yerr=mtBid_se, color='k', label='mtBid');
    errorbar(n_cpt_time,  n_cpt_means[dt['ctBid']] / n_cpt_tBid_0, \
             yerr=ctBid_se, color='m', label='ctBid');

    plot(time, one_cpt_obs['mBax'] / one_cpt_Bax_0, color='r');
    plot(time, one_cpt_obs['cBax'] / one_cpt_Bax_0, color='b');
    plot(time, one_cpt_obs['iBax'] / one_cpt_Bax_0, color='g');
    plot(time, one_cpt_obs['mtBid'] / one_cpt_tBid_0, color='k');
    plot(time, one_cpt_obs['ctBid'] / one_cpt_tBid_0, color='m');

    title('one_cpt vs. n_cpt')
    xlabel('Time (sec)')
    ylabel('Fraction of Total Bax/tBid')
    legend(loc='lower right')

if __name__ == '__main__':
    from __init__ import jobs
    f = h5py.File('data.hdf5')
    data = f['data']
    dtypes = pickle.loads(f['data_dtype_pck'][:])
    dtypes_dict = dict((name, i) for i, name in enumerate(dtypes.names))
    plot_timecourse_comparison(jobs[0], data[0,:,:,:], dtypes_dict)

