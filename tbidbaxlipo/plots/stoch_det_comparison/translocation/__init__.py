from tbidbaxlipo.models import simulation
from pylab import *

def plot_timecourse_comparison(job, n_cpt_obs):
    b_one = job.one_cpt_builder()
    [time, one_cpt_obs] = job.run_one_cpt()
    one_cpt_Bax_0 = b_one.model.parameters['Bax_0'].value

    b_n = job.n_cpt_builder()
    #n_cpt_obs = job.run_n_cpt(cleanup=False)
    [n_cpt_means, n_cpt_sd] = simulation.calculate_mean_and_std(n_cpt_obs)
    n_cpt_Bax_0 = b_n.model.parameters['Bax_0'].value

    # Compare deterministic to n_cpt
    figure()
    mBax_se = (n_cpt_sd['mBax'] / n_cpt_Bax_0) / np.sqrt(job.num_sims)
    cBax_se = (n_cpt_sd['cBax'] / n_cpt_Bax_0) / np.sqrt(job.num_sims)
    errorbar(n_cpt_means['time'], n_cpt_means['mBax'] / n_cpt_Bax_0, \
             yerr=mBax_se, color='r', label='mBax');
    errorbar(n_cpt_means['time'], n_cpt_means['cBax'] / n_cpt_Bax_0, \
             yerr=cBax_se, color='b', label='cBax');
    plot(time, one_cpt_obs['mBax'] / one_cpt_Bax_0, color='r');
    plot(time, one_cpt_obs['cBax'] / one_cpt_Bax_0, color='b');
    title('one_cpt vs. n_cpt')
    xlabel('Time (sec)')
    ylabel('Fraction of Total Bax')
    legend(loc='lower right')

    # Compare deterministic to site_cpt
    """
    b_site = job.site_cpt_builder()
    site_cpt_obs = job.run_site_cpt()
    [site_cpt_means, site_cpt_sd] = \
                        simulation.calculate_mean_and_std(site_cpt_obs)
    site_cpt_Bax_0 = b_site.model.parameters['Bax_0'].value
    figure()
    errorbar(site_cpt_means['time'], site_cpt_means['mBax'] / site_cpt_Bax_0, \
             yerr=site_cpt_sd['mBax'] / site_cpt_Bax_0, color='r', label='mBax')
    errorbar(site_cpt_means['time'], site_cpt_means['cBax'] / site_cpt_Bax_0, 
             yerr=site_cpt_sd['cBax'] / site_cpt_Bax_0, color='b', label='cBax')
    plot(time, one_cpt_obs['mBax'] / one_cpt_Bax_0, color='r');
    plot(time, one_cpt_obs['cBax'] / one_cpt_Bax_0, color='b');
    title('one_cpt vs. site_cpt')
    xlabel('Time (sec)')
    ylabel('Fraction of Total Bax')
    legend(loc='lower right')
    """

