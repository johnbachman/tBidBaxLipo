from tbidbaxlipo.models import simulation
import pickle
import pkgutil
import os
import sys

mod_path = os.path.dirname(sys.modules[__name__].__file__)
hdf5_filename = os.path.abspath(os.path.join(mod_path, 'data.hdf5'))
if os.path.exists(hdf5_filename):
    data = simulation.CptDataset(hdf5_filename)

class Job(simulation.Job):
    def __init__(self, tbid_conc):
        params_dict = {'tBid_transloc_kr': 0, 'tBid_0': tbid_conc}
        scaling_factor = 20
        tmax = 10000
        num_sims = 30
        n_steps = 100
        super(Job, self).__init__(params_dict, scaling_factor, tmax, n_steps,
                                  num_sims)

    def build(self, module):
        builder = module.Builder(params_dict=self.params_dict,
                                 scaling_factor=self.scaling_factor)
        builder.translocate_tBid_Bax()
        builder.tBid_activates_Bax()
        builder.pores_from_Bax_monomers()
        builder.model.name = job_name
        return builder

# Create jobs for each condition
tbid_concs = range(21)
jobs = [Job(tbid_conc) for tbid_conc in tbid_concs]
job_name = 'tbid_titration'

if __name__ == '__main__':
    simulation.run_main(jobs)
