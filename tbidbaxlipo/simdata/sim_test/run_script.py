from tbidbaxlipo.models import n_cpt_jobs

class Job(n_cpt_jobs.Job):
    def __init__(self):
        params_dict = {'Bax_transloc_kf':1e-2}
        scaling_factor = 5
        tmax = 10000
        n_steps = 200
        num_sims = 3
        super(Job, self).__init__(params_dict, scaling_factor, tmax,
                                  n_steps, num_sims)

    def build(self, module):
        builder = module.Builder(params_dict=self.params_dict,
                                 scaling_factor=self.scaling_factor)
        builder.build_model_bax_heat()
        return builder

if __name__ == '__main__':
    j = Job()
    j.run_n_cpt()
