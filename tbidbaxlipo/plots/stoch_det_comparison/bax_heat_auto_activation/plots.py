from pylab import *
import numpy as np

ion()

def plot_timecourse_comparison(jobs, data, cond_index):
    job = jobs[cond_index]
    b_one = job.one_cpt_builder()
    [time, one_cpt_obs] = job.run_one_cpt()
    one_cpt_Bax_0 = b_one.model.parameters['Bax_0'].value

    b_n = job.n_cpt_builder()
    n_cpt_Bax_0 = b_n.model.parameters['Bax_0'].value

    # Compare deterministic to n_cpt
    figure()
    mBax_se = ((data.sds(cond_index, 'mBax') / n_cpt_Bax_0) /
               np.sqrt(job.num_sims))
    cBax_se = ((data.sds(cond_index, 'cBax') / n_cpt_Bax_0) /
               np.sqrt(job.num_sims))
    iBax_se = ((data.sds(cond_index, 'iBax') / n_cpt_Bax_0) /
               np.sqrt(job.num_sims))
    iBax_se = ((data.sds(cond_index, 'pBax') / n_cpt_Bax_0) /
               np.sqrt(job.num_sims))

    cpt_time = data.means(cond_index, 'time')
    errorbar(cpt_time,  data.means(cond_index, 'mBax') / n_cpt_Bax_0, \
             yerr=mBax_se, color='r', label='mBax');
    errorbar(cpt_time,  data.means(cond_index, 'cBax') / n_cpt_Bax_0, \
             yerr=cBax_se, color='b', label='cBax');
    errorbar(cpt_time,  data.means(cond_index, 'iBax') / n_cpt_Bax_0, \
             yerr=iBax_se, color='g', label='iBax');
    errorbar(cpt_time,  data.means(cond_index, 'pBax') / n_cpt_Bax_0, \
             yerr=iBax_se, color='m', label='pBax');

    plot(time, one_cpt_obs['mBax'] / one_cpt_Bax_0, color='r');
    plot(time, one_cpt_obs['cBax'] / one_cpt_Bax_0, color='b');
    plot(time, one_cpt_obs['iBax'] / one_cpt_Bax_0, color='g');
    plot(time, one_cpt_obs['pBax'] / one_cpt_Bax_0, color='m');

    title('one_cpt vs. n_cpt')
    xlabel('Time (sec)')
    ylabel('Fraction of Total Bax')
    legend(loc='lower right')

if __name__ == '__main__':
    from tbidbaxlipo.plots.stoch_det_comparison.bax_heat_auto_activation \
            import jobs, data
    plot_timecourse_comparison(jobs, data, 0)
