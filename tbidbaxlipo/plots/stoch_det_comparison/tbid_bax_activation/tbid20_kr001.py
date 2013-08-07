from tbidbaxlipo.models import simulation
import pickle
import pkgutil

class Job(simulation.Job):
    def __init__(self):
        params_dict = {'tBid_transloc_kr': 1e-2, 'tBid_0':20}
        scaling_factor = 20
        tmax = 10000
        num_sims = 3
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

job = Job()
job_name = 'tbid20_kr001'

try:
    n_cpt_obs = pickle.loads(pkgutil.get_data(
                   'tbidbaxlipo.plots.stoch_det_comparison.tbid_bax_activation',
                   '%s.pck' % job_name))
except IOError:
    pass

if __name__ == '__main__':
    n_cpt_obs = job.run_n_cpt(cleanup=True)
    with open('%s.pck' % job_name, 'w') as f:
        pickle.dump(n_cpt_obs, f)

