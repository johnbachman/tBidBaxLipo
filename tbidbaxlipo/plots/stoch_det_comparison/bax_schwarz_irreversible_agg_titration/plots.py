from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

plt.ion()

def calc_min_pore_size(jobs, data):
    # For each concentration, get final dye release percentage, or final pore
    # num (start with dye release pctage)
    dr_fmax_means = np.zeros(len(jobs))
    dr_fmax_ses = np.zeros(len(jobs))
    bax_concs = np.zeros(len(jobs))
    num_sims = data.sim_data.shape[1]
    lipo_conc = jobs[0].one_cpt_builder().model.parameters['Vesicles_0'].value
    min_pore_sizes = np.zeros(len(jobs))

    for cond_index, job in enumerate(jobs):
        (dr_mean, dr_sd) = data.get_mean_dye_release(cond_index)
        dr_se = dr_sd / np.sqrt(num_sims)
        dr_fmax_means[cond_index] = dr_mean[-1]
        dr_fmax_ses[cond_index] = dr_se[-1]
        bax_concs[cond_index] = \
                       job.one_cpt_builder().model.parameters['Bax_0'].value
        ratio = bax_concs[cond_index] / float(lipo_conc)
        min_pore_sizes[cond_index] = \
                stats.poisson.isf(dr_mean[-1], ratio)
    plt.figure()
    plt.plot(bax_concs, min_pore_sizes, marker='o')
    plt.xlabel('[Bax]')
    plt.ylabel('Min pore size')

    lin_fit = stats.linregress(bax_concs, min_pore_sizes)
    slope = lin_fit[0]
    intercept = lin_fit[1]
    plt.plot(bax_concs, slope*bax_concs + intercept, color='r')

    plt.title('Slope %f, intercept %f' % (slope, intercept))
