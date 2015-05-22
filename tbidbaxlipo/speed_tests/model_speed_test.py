from tbidbaxlipo.models.nbd.multiconf import Builder
from pysb.integrate import Solver
import numpy as np
import timeit

bd = Builder()
bd.build_model_multiconf(2, 1., 0.1, 10, normalized_data=True)
bd.model.parameters['c1_scaling'].value = 5.

t = np.linspace(0, 4000, 150)
sol = Solver(bd.model, t)

num_simulations = 10000
start_time = timeit.default_timer()
for i in xrange(0, num_simulations):
    sol.run()
elapsed = timeit.default_timer() - start_time

print '%d simulations in %f seconds' % (num_simulations, elapsed)
print '%f per simulation' % (elapsed / float(num_simulations))

