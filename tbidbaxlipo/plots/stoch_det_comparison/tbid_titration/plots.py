from pylab import *

def plot_dye_release_titration(jobs, data):
    figure()
    num_conditions = len(jobs)
    cpt_time = data.sim_data[0,0,0,:]
    for cond_index in range(num_conditions):
        (dr_mean, dr_sd) = data.get_mean_dye_release(cond_index)
        errorbar(cpt_time, dr_mean, yerr=dr_sd, label='Cond %d' % cond_index)
    xlabel('Time (sec)')
    ylabel('Pct. Dye Release')
    title('Dye release kinetics for tBid titration')
    legend(loc='lower right')
    show()

if __name__ == '__main__':
    ion()
    from __init__ import jobs, data
    plot_dye_release_titration(jobs, data)
