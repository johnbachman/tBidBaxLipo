from tbidbaxlipo.models import simulation
import pickle
import os
import sys

mod_path = os.path.dirname(sys.modules[__name__].__file__)
hdf5_filename = os.path.abspath(os.path.join(mod_path, 'data.hdf5'))
if os.path.exists(hdf5_filename):
    data = simulation.CptDataset(hdf5_filename)

class Job(simulation.Job):
    def __init__(self, Bax_0):
        params_dict = {'Bax_0': Bax_0,
                       'tBid_0': 0,
                       'basal_Bax_kf': 1e-3, 'basal_Bax_kr': 1e-2,
                       'iBax_activates_mBax_k': 6e-3}
        scaling_factor = 5
        tmax = 10000
        num_sims = 3
        n_steps = 100
        super(Job, self).__init__(params_dict, scaling_factor, tmax, n_steps,
                                  num_sims)

    def build(self, module):
        builder = module.Builder(params_dict=self.params_dict,
                                 scaling_factor=self.scaling_factor)
        builder.build_model_bax_heat_auto_reversible_activation()
        return builder

# Create jobs for each condition
jobs = [Job(200)]
job_name = 'bax_heat_auto_titration'

if __name__ == '__main__':
    simulation.run_main(jobs)
