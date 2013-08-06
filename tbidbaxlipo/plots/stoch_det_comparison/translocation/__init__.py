from tbidbaxlipo.models import one_cpt, n_cpt, site_cpt, simulation
from pylab import *
import collections
import glob
from scipy.stats import poisson
import pkgutil
import pickle

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

ion()
job = Job()

n_cpt_obs = pickle.loads(pkgutil.get_data(
                    'tbidbaxlipo.plots.stoch_det_comparison.translocation',
                    'n_cpt_obs.pck'))

def plot_timecourse_comparison():
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

def plot_hist_vs_poisson(observable_basename, timepoint_index):
    # Get the things we'll need
    b_n = job.n_cpt_builder()
    b_one = job.one_cpt_builder()
    [time, one_cpt_obs] = job.run_one_cpt()

    # Get histogram data
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
    xlabel('%s per liposome' % observable_basename)
    ylabel('Frequency')

    return means

def print_mbax_means_and_vars(timepoint):
    b_n = job.n_cpt_builder()
    (mbax_means, mbax_vars) = simulation.get_means_across_cpts(
                        b_n.get_compartment_observables('mBax'),
                        n_cpt_obs, timepoint=timepoint)
    print "mBax per liposome mean: %f +/- %f" % \
          (np.mean(mbax_means), np.std(mbax_means))
    print "mBax per liposome variance: %f +/- %f" % \
          (np.mean(mbax_vars), np.std(mbax_vars))

if __name__ == '__main__':
    job.run_n_cpt(cleanup=False)
