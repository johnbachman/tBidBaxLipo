from tbidbaxlipo.models import simulation
from tbidbaxlipo.models import one_cpt
from matplotlib import pyplot as plt
import numpy as np

class Job(simulation.Job):
    def __init__(self):
        params_dict = {'Vesicles_0':50, 'Bax_0':50}
        scaling_factor = 1
        tmax = 10000
        n_steps = 200
        num_sims = 20
        super(Job, self).__init__(params_dict, scaling_factor, tmax,
                                  n_steps, num_sims)

    def build(self, module):
        builder = module.Builder(params_dict=self.params_dict,
                                 scaling_factor=self.scaling_factor)
        builder.build_model_bax_schwarz()
        return builder

def plot():
    # Create the job instance and get a builder
    j = Job()
    one_cpt_builder = j.build(one_cpt)

    # Run the deterministic simulation
    (t, det_obs) = j.run_one_cpt()

    # Run the stochastic simulation
    xrecs = j.run_n_cpt(cleanup=True)
    (means, stds) = simulation.calculate_mean_and_std(xrecs)
    (dr_mean, dr_std) = simulation.calculate_dye_release_mean_and_std(xrecs)

    plt.ion()
    plt.figure()

    # Pores
    plt.plot(t, det_obs['pores'], label='__nolabel__', color='r')
    plt.errorbar(means['time'], means['pores'] / j.scaling_factor,
                 yerr=stds['pores'] / j.scaling_factor, label='Pores',
                 color='r')
    # Bax_c
    plt.plot(t, det_obs['cBax'], color='g')
    plt.errorbar(means['time'], means['cBax'] / j.scaling_factor,
                 yerr=stds['pores'] / j.scaling_factor, label='cBax',
                 color='g')
    # Bax_m
    plt.plot(t, det_obs['mBax'], color='b')
    plt.errorbar(means['time'], means['mBax'] / j.scaling_factor,
                 yerr=stds['pores'] / j.scaling_factor, label='mBax',
                 color='b')

    # Label the plot
    plt.xlabel('Time (secs)')
    plt.ylabel('Total pores')
    plt.title('Comparing one_cpt and n_cpt simulations')
    plt.legend(loc='lower right')

    # -- Dye release --
    # Plot deterministic results
    plt.figure()
    plt.xlabel('Time (secs)')
    plt.ylabel('Dye release')
    num_vesicles = one_cpt_builder.model.parameters['Vesicles_0'].value
    avg_pores = det_obs['pores'] / num_vesicles
    plt.plot(t, 1 - np.exp(-avg_pores), color='b')
    plt.errorbar(means['time'], dr_mean, yerr=dr_std, color='b',
                 label='Dye release')
    plt.legend(loc='lower right')

if __name__ == '__main__':
    plot()
