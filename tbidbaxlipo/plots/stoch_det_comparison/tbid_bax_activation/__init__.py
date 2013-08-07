from tbidbaxlipo.models import one_cpt, n_cpt, site_cpt, simulation
from pylab import *
import collections
import glob
from scipy.stats import poisson
import pkgutil
import pickle

ion()

def plot_timecourse_comparison(job, n_cpt_obs):
    b_one = job.one_cpt_builder()
    [time, one_cpt_obs] = job.run_one_cpt()
    one_cpt_Bax_0 = b_one.model.parameters['Bax_0'].value
    one_cpt_tBid_0 = b_one.model.parameters['tBid_0'].value

    b_n = job.n_cpt_builder()
    [n_cpt_means, n_cpt_sd] = simulation.calculate_mean_and_std(n_cpt_obs)
    n_cpt_Bax_0 = b_n.model.parameters['Bax_0'].value
    n_cpt_tBid_0 = b_n.model.parameters['tBid_0'].value

    # Compare deterministic to n_cpt
    figure()
    mBax_se = (n_cpt_sd['mBax'] / n_cpt_Bax_0) / np.sqrt(job.num_sims)
    cBax_se = (n_cpt_sd['cBax'] / n_cpt_Bax_0) / np.sqrt(job.num_sims)
    iBax_se = (n_cpt_sd['iBax'] / n_cpt_Bax_0) / np.sqrt(job.num_sims)
    mtBid_se = (n_cpt_sd['mtBid'] / n_cpt_tBid_0) / np.sqrt(job.num_sims)
    ctBid_se = (n_cpt_sd['ctBid'] / n_cpt_tBid_0) / np.sqrt(job.num_sims)

    errorbar(n_cpt_means['time'], n_cpt_means['mBax'] / n_cpt_Bax_0, \
             yerr=mBax_se, color='r', label='mBax');
    errorbar(n_cpt_means['time'], n_cpt_means['cBax'] / n_cpt_Bax_0, \
             yerr=cBax_se, color='b', label='cBax');
    errorbar(n_cpt_means['time'], n_cpt_means['iBax'] / n_cpt_Bax_0, \
             yerr=iBax_se, color='g', label='iBax');
    errorbar(n_cpt_means['time'], n_cpt_means['mtBid'] / n_cpt_tBid_0, \
             yerr=mtBid_se, color='k', label='mtBid');
    errorbar(n_cpt_means['time'], n_cpt_means['ctBid'] / n_cpt_tBid_0, \
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

