from pylab import *
from tbidbaxlipo.models import simulation
from scipy.stats import poisson

def plot_bax_timecourse_comparison(jobs, data, cond_index):
    """Compare stochastic and deterministic timecourses for Bax observables.

    Parameters
    ----------
    jobs : list of Job objects
        Each job in the list represents a simulation condition and can be
        used to extract the simulation parameters.
    data : instance of CptDataset
        CptDataset object containing the simulation data drawn from the
        HDF5 datafile.
    cond_index : int
        The index of the condition for which we wish to plot the observable,
        i.e., the condition represented by the job in jobs[cond_index].
    """

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
    pBax_se = ((data.sds(cond_index, 'pBax') / n_cpt_Bax_0) /
               np.sqrt(job.num_sims))

    cpt_time = data.means(cond_index, 'time')
    errorbar(cpt_time,  data.means(cond_index, 'mBax') / n_cpt_Bax_0, \
             yerr=mBax_se, color='r', label='mBax');
    errorbar(cpt_time,  data.means(cond_index, 'cBax') / n_cpt_Bax_0, \
             yerr=cBax_se, color='b', label='cBax');
    errorbar(cpt_time,  data.means(cond_index, 'iBax') / n_cpt_Bax_0, \
             yerr=iBax_se, color='g', label='iBax');
    errorbar(cpt_time,  data.means(cond_index, 'pBax') / n_cpt_Bax_0, \
             yerr=pBax_se, color='m', label='pBax');

    plot(time, one_cpt_obs['mBax'] / one_cpt_Bax_0, color='r');
    plot(time, one_cpt_obs['cBax'] / one_cpt_Bax_0, color='b');
    plot(time, one_cpt_obs['iBax'] / one_cpt_Bax_0, color='g');
    plot(time, one_cpt_obs['pBax'] / one_cpt_Bax_0, color='m');

    title('one_cpt vs. n_cpt')
    xlabel('Time (sec)')
    ylabel('Fraction of Total Bax')
    legend(loc='lower right')

def plot_dye_release_titration(jobs, data):
    figure()
    num_sims = data.sim_data.shape[1]
    cpt_time = data.sim_data[0,0,0,:]
    for cond_index, job in enumerate(jobs):
        b_one = job.one_cpt_builder()
        (one_cpt_time, one_cpt_obs) = job.run_one_cpt()
        avg_pores = one_cpt_obs['pores'] / \
                    b_one.model.parameters['Vesicles_0'].value
        one_cpt_dr = 1 - np.exp(-avg_pores)
        (dr_mean, dr_sd) = data.get_mean_dye_release(cond_index)
        dr_se = dr_sd / sqrt(num_sims)
        plot(one_cpt_time, one_cpt_dr, color='b')
        errorbar(cpt_time, dr_mean, yerr=dr_se, color='r',
                 label='Cond %d' % cond_index)
    xlabel('Time (sec)')
    ylabel('Pct. Dye Release')
    title('Dye release kinetics')
    #legend(loc='lower right')
    show()

def plot_dr_fmax_vs_bax(jobs, data):
    dr_fmax_means = np.zeros(len(jobs))
    dr_fmax_ses = np.zeros(len(jobs))
    bax_concs = np.zeros(len(jobs))
    num_sims = data.sim_data.shape[1]
    for cond_index, job in enumerate(jobs):
        (dr_mean, dr_sd) = data.get_mean_dye_release(cond_index)
        dr_se = dr_sd / sqrt(num_sims)
        dr_fmax_means[cond_index] = dr_mean[-1]
        dr_fmax_ses[cond_index] = dr_se[-1]
        bax_concs[cond_index] = \
                       job.one_cpt_builder().model.parameters['Bax_0'].value

    plt.figure()
    plt.plot(bax_concs, dr_fmax_means)
    plt.xlabel('[Bax]')
    plt.ylabel('Final dye release')
    plt.title('Max dye release vs. [Bax]')
    plt.show()

    plt.figure()
    plt.plot(bax_concs, -np.log10(1 - dr_fmax_means))
    plt.xlabel('[Bax]')
    plt.ylabel('Inferred max pores')
    plt.title('Max pores vs. [Bax]')
    plt.show()

def plot_hist_vs_poisson(jobs, data, cond_index, obs_basename, timepoint_index):
    """Plot the compartment distribution of an observable against a Poisson.

    Plots the distribution of the number of compartments that have a certain
    amount of the observable of interest, and compares this distribution to a
    Poisson distribution with the mean drawn from the comparable deterministic
    simulation.

    To make the plot, the frequency distribution for each simulation is
    computed, and then the mean frequency at each value is take across the
    simulations. The results are presented as a bar plot with error bars
    representing the standard error of the frequency over the simulations.

    Parameters
    ----------
    jobs : list of Job objects
        Each job in the list represents a simulation condition and can be
        used to extract the simulation parameters.
    data : instance of CptDataset
        CptDataset object containing the simulation data drawn from the
        HDF5 datafile.
    cond_index : int
        The index of the condition for which we wish to plot the observable,
        i.e., the condition represented by the job in jobs[cond_index].
    obs_basename : string
        The base name of the observable of interest. Observables for each
        compartment are expected to follow the pattern of obs_basename +
        '_c1', '_c2', '_c3'. Note that the obs_basename should not include
        the underscore separate the basename from the compartment identifier.
    timepoint_index : int
        The timepoint at which we wish to plot the distributions. Ranges
        from 0 to jobs[cond_index].n_steps.
    """

    job = jobs[cond_index]
    # Get the things we'll need
    b_n = job.n_cpt_builder()
    b_one = job.one_cpt_builder()
    [time, one_cpt_obs] = job.run_one_cpt()

    # Get histogram data
    (index, freq_matrix) = data.get_frequency_matrix(cond_index,
                                obs_basename, timepoint_index)
    means = np.mean(freq_matrix, axis=1)
    sds = np.std(freq_matrix, axis=1)

    norm_means = means / np.sum(means)
    norm_sds = sds / np.sum(means)
    norm_se = norm_sds / np.sqrt(freq_matrix.shape[1])

    figure()
    plot_index = index - 0.5 # offset so middle of bars line up with ints
    bar(plot_index, norm_means, 1.0, yerr=norm_se, color='r')

    # Value in deterministic model
    obs_per_lipo = one_cpt_obs[obs_basename][timepoint_index] / \
                    b_one.model.parameters['Vesicles_0'].value

    plot(index, poisson.pmf(index, obs_per_lipo), color='g',
         linewidth='2')
    xlim([min(index)-0.5, max(index)+0.5])
    title('Distribution of %s at %.1f seconds' %
          (obs_basename, time[timepoint_index]))
    xlabel('%s per liposome' % obs_basename)
    ylabel('Frequency')

    return means

def print_obs_means_and_vars(jobs, data, cond_index, obs_basename, timepoint):
    job = jobs[cond_index]
    b_n = job.n_cpt_builder()
    (obs_means, obs_vars) = data.get_means_across_cpts(cond_index,
                                                       obs_basename, timepoint)
    print "%s per liposome mean: %f +/- %f" % \
          (obs_basename, np.mean(obs_means), np.std(obs_means))
    print "%s per liposome variance: %f +/- %f" % \
          (obs_basename, np.mean(obs_vars), np.std(obs_vars))

