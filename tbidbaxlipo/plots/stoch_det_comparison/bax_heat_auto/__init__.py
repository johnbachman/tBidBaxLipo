from tbidbaxlipo.models import simulation
import pickle
import os
import sys

mod_path = os.path.dirname(sys.modules[__name__].__file__)
hdf5_filename = os.path.abspath(os.path.join(mod_path, 'data.hdf5'))
if os.path.exists(hdf5_filename):
    data = simulation.CptDataset(hdf5_filename)

class Job(simulation.Job):
    def __init__(self):
        params_dict = {'tBid_transloc_kr': 0, 'tBid_0': 0}
        scaling_factor = 20
        tmax = 10000
        num_sims = 30
        n_steps = 100
        super(Job, self).__init__(params_dict, scaling_factor, tmax, n_steps,
                                  num_sims)

    def build(self, module):
        builder = module.Builder(params_dict=self.params_dict,
                                 scaling_factor=self.scaling_factor)
        builder.build_model_bax_heat_auto1()
        return builder

# Create jobs for each condition
jobs = [Job()]
job_name = 'bax_heat_auto1'

if __name__ == '__main__':
    simulation.run_main(jobs)
