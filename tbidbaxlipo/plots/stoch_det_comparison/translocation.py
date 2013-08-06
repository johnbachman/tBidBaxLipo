from tbidbaxlipo.models import one_cpt, n_cpt, site_cpt, simulation
from pylab import *
import collections

class Job(simulation.Job):
    def __init__(self):
        scaling_factor = 10
        tmax = 60
        num_sims = 40
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
    n_cpt_obs = job.run_n_cpt(cleanup=True)
    [n_cpt_means, n_cpt_sd] = job.calculate_mean_and_std(n_cpt_obs)
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
    [site_cpt_means, site_cpt_sd] = job.calculate_mean_and_std(site_cpt_obs)
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

    # Look at histograms
    counts_list = []
    for obs in n_cpt_obs:
        values = job.get_observables_values(
                        b_n.get_compartment_observables('mBax'), obs)
        counts_list.append(collections.Counter(values))
        #figure()
        #bins = range(int(min(values)), int(max(values)))
        #hist(counts, bins=bins)

    all_keys = set()
    for counts in counts_list:
        all_keys |= set(counts.keys())

    key_min = min(all_keys)
    key_max = max(all_keys)
    counts_arr = np.zeros((key_max - key_min + 1, len(counts_list)))

    for i, counts in enumerate(counts_list):
        for key, val in counts.iteritems():
            counts_arr[key - key_min, i] = val

    means = np.mean(counts_arr, axis=1)
    sds = np.std(counts_arr, axis=1)
    ind = range(int(key_min), int(key_max) + 1)
    figure()
    bar(ind, means, 0.6, yerr=sds, color='r')
