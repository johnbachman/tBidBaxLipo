from pylab import *

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
    # Plot n_cpt model
    cpt_time = data.means(cond_index, 'time')
    errorbar(cpt_time,  data.means(cond_index, 'mBax') / n_cpt_Bax_0, \
             yerr=mBax_se, color='r', label='mBax');
    errorbar(cpt_time,  data.means(cond_index, 'cBax') / n_cpt_Bax_0, \
             yerr=cBax_se, color='b', label='cBax');
    # Plot deterministic model
    plot(time, one_cpt_obs['mBax'] / one_cpt_Bax_0, color='r');
    plot(time, one_cpt_obs['cBax'] / one_cpt_Bax_0, color='b');

    title('one_cpt vs. n_cpt')
    xlabel('Time (sec)')
    ylabel('Fraction of Total Bax')
    legend(loc='lower right')

if __name__ == '__main__':
    from __init__ import jobs, data
    plot_timecourse_comparison(jobs, data, 0)
