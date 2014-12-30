from pysb.integrate import Solver
import numpy as np

from tbidbaxlipo.models.one_cpt import Builder

bd = Builder()

bd.build_model_nbd_2_conf_dimer()

t = np.linspace(0, 1e4, 1e3)
sol = Solver(bd.model, t)
sol.run()
