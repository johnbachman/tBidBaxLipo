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
        params_dict = {'Bax_0': Bax_0,
                       'tBid_0': 0.,
                       'Vesicles_0': 5.,
                       'Bax_transloc_kf': 1e-3,
                       'Bax_transloc_kr': 1e-1,
                       'pore_formation_rate_k': 1e-3,
                       }
        scaling_factor = 20
        tmax = 1000
        num_sims = 30
        n_steps = 100
        super(Job, self).__init__(params_dict, scaling_factor, tmax, n_steps,
                                  num_sims)

    def build(self, module):
        builder = module.Builder(params_dict=self.params_dict,
                                 scaling_factor=self.scaling_factor)
        builder.build_model_bax_schwarz_irreversible()
        return builder

# Create jobs for each condition
# Make sure concs are integers for best agreement between 1C and MC models
bax_concs = np.round(np.logspace(0, 3, 20))
jobs = [Job(bax_conc) for bax_conc in bax_concs]
job_name = 'bax_schwarz_irreversible_titration'

if __name__ == '__main__':
    simulation.run_main(jobs)
