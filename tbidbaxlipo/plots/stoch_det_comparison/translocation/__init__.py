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
        params_dict = {}
        scaling_factor = 10
        tmax = 60
        num_sims = 60
        n_steps = 100
        super(Job, self).__init__(params_dict, scaling_factor, tmax, n_steps,
                                  num_sims)

    def build(self, module):
        builder = module.Builder(params_dict=self.params_dict,
                                 scaling_factor=self.scaling_factor)
        builder.translocate_Bax()
        builder.model.name = job_name
        return builder

jobs = [Job()]
job_name = 'translocation'

if __name__ == '__main__':
    simulation.run_main(jobs)
