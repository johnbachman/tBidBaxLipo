from pylab import *

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
        plot(one_cpt_time, one_cpt_dr)
        errorbar(cpt_time, dr_mean, yerr=dr_se, label='Cond %d' % cond_index)
        xlim([0, 4000])
    xlabel('Time (sec)')
    ylabel('Pct. Dye Release')
    title('Dye release kinetics for tBid titration')
    legend(loc='lower right')
    show()

if __name__ == '__main__':
    ion()
    from __init__ import jobs, data
    plot_dye_release_titration(jobs, data)
