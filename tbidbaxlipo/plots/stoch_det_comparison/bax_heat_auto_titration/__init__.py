from tbidbaxlipo.models import simulation
import pickle
import os
import sys
from tbidbaxlipo.plots import layout_130614
import numpy as np

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
        scaling_factor = 20
        tmax = 10000
        num_sims = 30
        n_steps = 100
        super(Job, self).__init__(params_dict, scaling_factor, tmax, n_steps,
                                  num_sims)

    def build(self, module):
        builder = module.Builder(params_dict=self.params_dict,
                                 scaling_factor=self.scaling_factor)
        builder.build_model_bax_heat_auto1_reversible_activation()
        return builder

# Create jobs for each condition
bax_concs = np.array(layout_130614.df.columns, dtype='float')
jobs = [Job(bax_conc) for bax_conc in bax_concs]
job_name = 'bax_heat_auto1_reversible_activation_titration'

if __name__ == '__main__':
    simulation.run_main(jobs)
