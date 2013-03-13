from pysb.integrate import Solver
import pymc
from nbd_parallel_model import model as nbd_model
import numpy as np
import nbd_analysis as nbd

tspan = nbd.time_other
nbd_avgs, nbd_stds = nbd.calc_norm_avg_std()
# Cut off the last two timepoints from C62 so that it's the same length
nbd_avgs[1] = nbd_avgs[1][0:-2]
nbd_stds[1] = nbd_stds[1][0:-2]

# Reduced data
#nbd_c3_reduced_avgs = nbd_avgs[0][::100]
#nbd_c3_reduced_stds = nbd_stds[0][::100]
#tspan_reduced = tspan[::100]

# -- Set the prior distributions --
# For every scaling parameter, initialize a normal distribution
c3_scaling = pymc.Normal('c3_scaling', mu=1., tau=(1./0.1), value=0.8)

# For every rate parameter, initialize a lognormal distribution
c3_insertion_rate = pymc.Lognormal('c3_insertion_rate', mu=np.log(1e-3),
                                tau=(1./(np.log(10)*np.log(1e1))), value=1e-2)

# -- Define the model --
solver = Solver(nbd_model, tspan)
#solver = Solver(nbd_model, tspan_reduced)
params = np.array([p.value for p in nbd_model.parameters])

@pymc.stochastic(plot=False, observed=True)
def nbd_pymc_model(c3_scaling=c3_scaling, c3_insertion_rate=c3_insertion_rate,
                   value=0):
    # Update param array and run
    params[1] = c3_scaling
    params[6] = c3_insertion_rate
    solver.run(params)
    c3_model = c3_scaling * solver.yobs['Baxc3']
    return -np.sum(((nbd_avgs[0] - c3_model)**2) / (2 * nbd_stds[0]**2))

    # Return just the c3 trajectory for now
    #return c3_model

# -- Define the data as a multivariate normal variable --
#tau = np.eye(len(tspan)) * (1 / nbd_stds[0]**2)
#tau = np.eye(len(tspan_reduced)) * (1 / nbd_stds[0]**2)

#output = pymc.MvNormal('output', mu=nbd_pymc_model, tau=tau, observed=True,
#                  value=nbd_c3_reduced_avgs)

