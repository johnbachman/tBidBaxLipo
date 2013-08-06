from tbidbaxlipo.models import simulation
import pickle
import pkgutil

class Job(simulation.Job):
    def __init__(self):
        params_dict = {'Bax_transloc_kr':1e-1}
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

job = Job()
job_name = 'bax_kr_nominal'

try:
    n_cpt_obs = pickle.loads(pkgutil.get_data(
                        'tbidbaxlipo.plots.stoch_det_comparison.translocation',
                        '%s.pck' % job_name))
except IOError:
    pass

if __name__ == '__main__':
    n_cpt_obs = job.run_n_cpt(cleanup=True)
    with open('%s.pck' % job_name, 'w') as f:
        pickle.dump(n_cpt_obs, f)

