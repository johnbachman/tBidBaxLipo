from tbidbaxlipo.models import simulation
import pickle
import sys
import os
import h5py

mod_path = os.path.dirname(sys.modules[__name__].__file__)
hdf5_filename = os.path.abspath(os.path.join(mod_path, 'data.hdf5'))
if os.path.exists(hdf5_filename):
    data = simulation.CptDataset(hdf5_filename)

class Job(simulation.Job):
    def __init__(self, tBid_0, Bax_0):
        params_dict = {'tBid_0': tBid_0,
                       'Bax_0': Bax_0,
                       'Vesicles_0': 1.,
                       'tBid_transloc_kf': 1e-3,
                       'tBid_transloc_kr': 1e-2,
                       'Bax_transloc_kf': 1e-3,
                       'Bax_transloc_kr': 1e-3,
                       'tBid_Bax_mem_bh3_kf': 1e-2,
                       'tBid_Bax_mem_bh3_kr': 1e-2,
                       }
        scaling_factor = 20
        tmax = 2000
        num_sims = 30 # FIXME
        n_steps = 100
        super(Job, self).__init__(params_dict, scaling_factor, tmax, n_steps,
                                  num_sims)

    def build(self, module):
        builder = module.Builder(params_dict=self.params_dict,
                                 scaling_factor=self.scaling_factor)
        #builder.translocate_tBid_Bax()
        #builder.tBid_binds_Bax()
        builder.build_model_bid_bax_fret()
        builder.model.name = job_name
        return builder

jobs = [Job(20., 100.)]
job_name = 'tbid_bax_binding'

if __name__ == '__main__':
    simulation.run_main(jobs)
