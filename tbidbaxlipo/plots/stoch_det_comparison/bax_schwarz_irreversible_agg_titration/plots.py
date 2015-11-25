from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis, \
    fontsize

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
    ratios = []
    for cond_index, job in enumerate(jobs):
        (dr_mean, dr_sd) = data.get_mean_dye_release(cond_index)
        dr_se = dr_sd / np.sqrt(num_sims)
        dr_fmax_means[cond_index] = dr_mean[-1]
        dr_fmax_ses[cond_index] = dr_se[-1]
        bax_concs[cond_index] = \
                       job.one_cpt_builder().model.parameters['Bax_0'].value
        ratio = bax_concs[cond_index] / float(lipo_conc)
        ratios.append(ratio)
        min_pore_sizes[cond_index] = \
                stats.poisson.isf(dr_mean[-1], ratio)
    ratios = np.array(ratios)
    fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
    ax = fig.gca()
    ax.plot(ratios, min_pore_sizes, marker='o', linestyle='', color='b',
            markersize=2)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('[Bax]/[Lipo]')
    ax.set_ylabel('Predicted pore size')
    ax.set_ylim(0.7, 250)
    ax.set_xlim(0.1, 300)
    lin_fit = stats.linregress(ratios[0:], min_pore_sizes[0:])
    slope = lin_fit[0]
    intercept = lin_fit[1]
    lbound = 0.1
    ubound = 200
    ax.plot(ratios, slope*ratios + intercept, color='r')
    ax.set_title('Slope %f, intercept %f' % (slope, intercept),
                 fontsize=fontsize)
    format_axis(ax)
    import ipdb; ipdb.set_trace()

if __name__ == '__main__':
    from tbidbaxlipo.plots.stoch_det_comparison.bax_schwarz_irreversible_agg_titration import jobs, data

    calc_min_pore_size(jobs, data)
