
from ssa_job import Job
from tbidbaxlipo.models.one_cpt import Builder

b = Builder()
b.build_model_bax_heat()

j = Job(b.model, 10000, 200, 10, 1)
j.run()
