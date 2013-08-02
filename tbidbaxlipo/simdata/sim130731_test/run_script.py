from tbidbaxlipo.models.n_cpt_jobs import Job
from tbidbaxlipo.models.n_cpt import Builder

scaling_factor = 10
b = Builder(scaling_factor=scaling_factor)

b.build_model_bax_heat()

tmax = 10000
n_steps = 200
num_sims = 10
j = Job(b.model, tmax, n_steps, num_sims)
j.run()
