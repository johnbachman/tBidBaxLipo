from tbidbaxlipo.models import simulation
import pickle
import os
import sys
import numpy as np

mod_path = os.path.dirname(sys.modules[__name__].__file__)
hdf5_filename = os.path.abspath(os.path.join(mod_path, 'data.hdf5'))
if os.path.exists(hdf5_filename):
    data = simulation.CptDataset(hdf5_filename)

class Job(simulation.Job):
    def __init__(self, Bax_0):
        params_dict = {'Bax_0':Bax_0,
                       'Bax_transloc_kr': 5e-2,
                       'Bax_transloc_kf': 1e-3,
                       'Vesicles_0': 5
                       }
        scaling_factor = 10
        tmax = 100
        num_sims = 30
        n_steps = 100
        super(Job, self).__init__(params_dict, scaling_factor, tmax, n_steps,
                                  num_sims)

    def build(self, module):
        builder = module.Builder(params_dict=self.params_dict,
                                 scaling_factor=self.scaling_factor)
        builder.translocate_Bax()
        builder.model.name = job_name
        return builder

bax_concs = np.linspace(1, 600, 30)
jobs = [Job(bax_conc) for bax_conc in bax_concs]
job_name = 'translocation_titration'

if __name__ == '__main__':
    simulation.run_main(jobs)
