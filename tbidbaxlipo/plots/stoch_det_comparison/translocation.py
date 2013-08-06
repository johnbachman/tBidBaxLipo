from tbidbaxlipo.models import one_cpt, n_cpt, site_cpt, simulation
from pylab import *
import collections
import glob
from scipy.stats import poisson

class Job(simulation.Job):
    def __init__(self):
        scaling_factor = 10
        tmax = 60
        num_sims = 60
        n_steps = 100
        super(Job, self).__init__({}, scaling_factor, tmax, n_steps, num_sims)

    def build(self, module):
        builder = module.Builder(params_dict=self.params_dict,
                                 scaling_factor=self.scaling_factor)
        builder.translocate_Bax()
        return builder

if __name__ == '__main__':
    job = Job()
    ion()

    b_one = job.build(one_cpt)
    [time, one_cpt_obs] = job.run_one_cpt()
    one_cpt_Bax_0 = b_one.model.parameters['Bax_0'].value

    b_n = job.build(n_cpt)
    n_cpt_obs = simulation.load_bng_files(glob.glob('simdata/*.gdat'))
    #n_cpt_obs = job.run_n_cpt(cleanup=False)
    [n_cpt_means, n_cpt_sd] = simulation.calculate_mean_and_std(n_cpt_obs)
    n_cpt_Bax_0 = b_n.model.parameters['Bax_0'].value

    # Compare deterministic to n_cpt
    figure()
    errorbar(n_cpt_means['time'], n_cpt_means['mBax'] / n_cpt_Bax_0, \
             yerr=n_cpt_sd['mBax'] / n_cpt_Bax_0, color='r', label='mBax');
    errorbar(n_cpt_means['time'], n_cpt_means['cBax'] / n_cpt_Bax_0, \
             yerr=n_cpt_sd['cBax'] / n_cpt_Bax_0, color='b', label='cBax');
    plot(time, one_cpt_obs['mBax'] / one_cpt_Bax_0, color='r');
    plot(time, one_cpt_obs['cBax'] / one_cpt_Bax_0, color='b');
    title('one_cpt vs. n_cpt')
    xlabel('Time (sec)')
    ylabel('Fraction of Total Bax')
    legend(loc='lower right')

    # Compare deterministic to site_cpt
    """
    b_site = job.build(site_cpt)
    site_cpt_obs = job.run_site_cpt()
    [site_cpt_means, site_cpt_sd] = \
                        simulation.calculate_mean_and_std(site_cpt_obs)
    site_cpt_Bax_0 = b_site.model.parameters['Bax_0'].value
    figure()
    errorbar(site_cpt_means['time'], site_cpt_means['mBax'] / site_cpt_Bax_0, \
             yerr=site_cpt_sd['mBax'] / site_cpt_Bax_0, color='r', label='mBax');
    errorbar(site_cpt_means['time'], site_cpt_means['cBax'] / site_cpt_Bax_0, \
             yerr=site_cpt_sd['cBax'] / site_cpt_Bax_0, color='b', label='cBax');
    plot(time, one_cpt_obs['mBax'] / one_cpt_Bax_0, color='r');
    plot(time, one_cpt_obs['cBax'] / one_cpt_Bax_0, color='b');
    title('one_cpt vs. site_cpt')
    xlabel('Time (sec)')
    ylabel('Fraction of Total Bax')
    legend(loc='lower right')
    """

    def plot_hist_vs_poisson(observable_basename, timepoint_index):
        (index, freq_matrix) = simulation.get_frequency_matrix(
                                    b_n.get_compartment_observables('mBax'),
                                    n_cpt_obs, timepoint=timepoint_index)
        means = np.mean(freq_matrix, axis=1)
        sds = np.std(freq_matrix, axis=1)

        norm_means = means / trapz(means, index)
        norm_sds = sds / trapz(means, index)
        norm_se = norm_sds / np.sqrt(len(n_cpt_obs))

        figure()
        plot_index = index - 0.5 # offset so middle of bars line up with ints
        bar(plot_index, norm_means, 1.0, yerr=norm_se, color='r')

        # Value in deterministic model
        mBax_per_lipo = one_cpt_obs['mBax'][timepoint_index] / \
                        b_one.model.parameters['Vesicles_0'].value

        plot(index, poisson.pmf(index, mBax_per_lipo), color='g',
             linewidth='2')
        xlim([min(plot_index), max(index)])
        title('Distribution of %s at %.1f seconds' %
              (observable_basename, time[timepoint_index]))

    plot_hist_vs_poisson('mBax', 1)
    plot_hist_vs_poisson('mBax', 10)
    plot_hist_vs_poisson('mBax', 20)
    plot_hist_vs_poisson('mBax', 30)
    plot_hist_vs_poisson('mBax', 40)
    plot_hist_vs_poisson('mBax', job.n_steps)


