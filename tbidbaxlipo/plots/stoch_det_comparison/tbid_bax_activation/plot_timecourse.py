from tbidbaxlipo.models import one_cpt, n_cpt, site_cpt, simulation
from pylab import *
import collections
import glob
from scipy.stats import poisson
import pkgutil
import pickle
import numpy as np

ion()

def plot_timecourse_comparison(job, data):
    b_one = job.one_cpt_builder()
    [time, one_cpt_obs] = job.run_one_cpt()
    one_cpt_Bax_0 = b_one.model.parameters['Bax_0'].value
    one_cpt_tBid_0 = b_one.model.parameters['tBid_0'].value

    b_n = job.n_cpt_builder()
    n_cpt_Bax_0 = b_n.model.parameters['Bax_0'].value
    n_cpt_tBid_0 = b_n.model.parameters['tBid_0'].value

    # Compare deterministic to n_cpt
    figure()
    mBax_se = (data.sds(0, 'mBax') / n_cpt_Bax_0) / np.sqrt(job.num_sims)
    cBax_se = (data.sds(0, 'cBax') / n_cpt_Bax_0) / np.sqrt(job.num_sims)
    iBax_se = (data.sds(0, 'iBax') / n_cpt_Bax_0) / np.sqrt(job.num_sims)
    ctBid_se = (data.sds(0, 'ctBid') / n_cpt_tBid_0) / np.sqrt(job.num_sims)
    mtBid_se = (data.sds(0, 'mtBid') / n_cpt_tBid_0) / np.sqrt(job.num_sims)

    cpt_time = data.means(0, 'time')
    errorbar(cpt_time,  data.means(0, 'mBax') / n_cpt_Bax_0, \
             yerr=mBax_se, color='r', label='mBax');
    errorbar(cpt_time,  data.means(0, 'cBax') / n_cpt_Bax_0, \
             yerr=cBax_se, color='b', label='cBax');
    errorbar(cpt_time,  data.means(0, 'iBax') / n_cpt_Bax_0, \
             yerr=iBax_se, color='g', label='iBax');
    errorbar(cpt_time,  data.means(0, 'mtBid') / n_cpt_tBid_0, \
             yerr=mtBid_se, color='k', label='mtBid');
    errorbar(cpt_time,  data.means(0, 'ctBid') / n_cpt_tBid_0, \
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
    from __init__ import jobs, data
    plot_timecourse_comparison(jobs[0], data)

