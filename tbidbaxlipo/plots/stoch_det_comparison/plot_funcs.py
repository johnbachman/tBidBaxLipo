from pylab import *
from tbidbaxlipo.models import simulation
from scipy.stats import poisson

def plot_hist_vs_poisson(job, n_cpt_obs, observable_basename, timepoint_index):
    # Get the things we'll need
    b_n = job.n_cpt_builder()
    b_one = job.one_cpt_builder()
    [time, one_cpt_obs] = job.run_one_cpt()

    # Get histogram data
    (index, freq_matrix) = simulation.get_frequency_matrix(
                           b_n.get_compartment_observables(observable_basename),
                           n_cpt_obs, timepoint=timepoint_index)
    means = np.mean(freq_matrix, axis=1)
    sds = np.std(freq_matrix, axis=1)

    norm_means = means / np.sum(means)
    norm_sds = sds / np.sum(means)
    norm_se = norm_sds / np.sqrt(len(n_cpt_obs))

    figure()
    plot_index = index - 0.5 # offset so middle of bars line up with ints
    bar(plot_index, norm_means, 1.0, yerr=norm_se, color='r')

    # Value in deterministic model
    obs_per_lipo = one_cpt_obs[observable_basename][timepoint_index] / \
                    b_one.model.parameters['Vesicles_0'].value

    plot(index, poisson.pmf(index, obs_per_lipo), color='g',
         linewidth='2')
    xlim([min(index)-0.5, max(index)+0.5])
    title('Distribution of %s at %.1f seconds' %
          (observable_basename, time[timepoint_index]))
    xlabel('%s per liposome' % observable_basename)
    ylabel('Frequency')

    return means

def print_obs_means_and_vars(job, n_cpt_obs, observable_basename, timepoint):
    b_n = job.n_cpt_builder()
    (obs_means, obs_vars) = simulation.get_means_across_cpts(
                        b_n.get_compartment_observables(observable_basename),
                        n_cpt_obs, timepoint=timepoint)
    print "%s per liposome mean: %f +/- %f" % \
          (observable_basename, np.mean(obs_means), np.std(obs_means))
    print "%s per liposome variance: %f +/- %f" % \
          (observable_basename, np.mean(obs_vars), np.std(obs_vars))

