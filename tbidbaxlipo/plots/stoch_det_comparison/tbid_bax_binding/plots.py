from pylab import *
import numpy as np

ion()

def plot_timecourse_comparison(jobs, data, cond_index):
    job = jobs[cond_index]
    b_one = job.one_cpt_builder()
    [time, one_cpt_obs] = job.run_one_cpt()
    one_cpt_Bax_0 = b_one.model.parameters['Bax_0'].value
    one_cpt_tBid_0 = b_one.model.parameters['tBid_0'].value

    b_n = job.n_cpt_builder()
    n_cpt_Bax_0 = b_n.model.parameters['Bax_0'].value
    n_cpt_tBid_0 = b_n.model.parameters['tBid_0'].value

    # Compare deterministic to n_cpt
    figure()
    mBax_se = ((data.sds(cond_index, 'mBax') / n_cpt_Bax_0) /
               np.sqrt(job.num_sims))
    cBax_se = ((data.sds(cond_index, 'cBax') / n_cpt_Bax_0) /
               np.sqrt(job.num_sims))
    ctBid_se = ((data.sds(cond_index, 'ctBid') / n_cpt_tBid_0) /
               np.sqrt(job.num_sims))
    mtBid_se = ((data.sds(cond_index, 'mtBid') / n_cpt_tBid_0) /
               np.sqrt(job.num_sims))
    frac_tBid_se = ((data.sds(cond_index, 'tBidBax') / n_cpt_tBid_0) /
               np.sqrt(job.num_sims))
    frac_Bax_se = ((data.sds(cond_index, 'tBidBax') / n_cpt_Bax_0) /
               np.sqrt(job.num_sims))

    cpt_time = data.means(cond_index, 'time')
    errorbar(cpt_time,  data.means(cond_index, 'cBax') / n_cpt_Bax_0, \
             yerr=cBax_se, color='r', label='cBax', linestyle='');
    errorbar(cpt_time,  data.means(cond_index, 'mBax') / n_cpt_Bax_0, \
             yerr=mBax_se, color='g', label='mBax', linestyle='');
    errorbar(cpt_time,  data.means(cond_index, 'tBidBax') / n_cpt_Bax_0, \
             yerr=frac_Bax_se, color='b', label='Bax bound', linestyle='');
    errorbar(cpt_time,  data.means(cond_index, 'ctBid') / n_cpt_tBid_0, \
             yerr=ctBid_se, color='r', label='ctBid', linestyle='');
    errorbar(cpt_time,  data.means(cond_index, 'mtBid') / n_cpt_tBid_0, \
             yerr=mtBid_se, color='g', label='mtBid', linestyle='');
    errorbar(cpt_time,  data.means(cond_index, 'tBidBax') / n_cpt_tBid_0, \
             yerr=frac_tBid_se, color='b', label='tBid bound', linestyle='');

    plot(time, one_cpt_obs['cBax'] / one_cpt_Bax_0, color='r');
    plot(time, one_cpt_obs['mBax'] / one_cpt_Bax_0, color='g');
    plot(time, one_cpt_obs['tBidBax'] / one_cpt_Bax_0, color='b');
    plot(time, one_cpt_obs['ctBid'] / one_cpt_tBid_0, color='r',
            linestyle='--');
    plot(time, one_cpt_obs['mtBid'] / one_cpt_tBid_0, color='g',
            linestyle='--');
    plot(time, one_cpt_obs['tBidBax'] / one_cpt_tBid_0, color='b',
            linestyle='--');

    title('one_cpt vs. n_cpt')
    xlabel('Time (sec)')
    ylabel('Fraction of Total Bax/tBid')
    legend(loc='lower right')
