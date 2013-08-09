from tbidbaxlipo.models import simulation
import pickle
import pkgutil
import sys
import os

class Job(simulation.Job):
    def __init__(self, tBid_0, tBid_transloc_kr):
        params_dict = {'tBid_0': tBid_0, 'tBid_transloc_kr': tBid_transloc_kr}
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
        builder.model.name = job_name
        return builder

jobs = [Job(1, 0), Job(1, 1e-2), Job(20, 0), Job(20, 1e-2)]
job_name = 'tbid_bax_activation'

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "The condition index must also be specified."
        sys.exit()
    job_index = int(sys.argv[1])
    job = jobs[job_index]
    data_dir = os.path.join(os.getcwd(), 'data_%d' % job_index)
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    os.chdir(data_dir)
    job.run_n_cpt(cleanup=False)