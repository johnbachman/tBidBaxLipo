"""A worked example of using emcee to fit a simple PySB model by MCMC. """

from pysb import *
import numpy as np
from pysb.integrate import odesolve, Solver
from matplotlib import pyplot as plt
import emcee
from matplotlib import pyplot as plt

# Define the model: simple exponential decay from 1 to 0
Model()
Monomer('A')
Parameter('A_0', 1)
Parameter('k', 0.2)
Rule('A_decay', A() >> None, k)
Initial(A(), A_0)
Observable('A_obs', A())

# Simulate the model to generate the synthetic data
ntimes = 100
tspan = np.linspace(0, 40, ntimes)
ysim = odesolve(model, tspan)

# Add error to the underlying data by adding a Gaussian value with a standard
# deviation of 0.1
seed = 1
rs = np.random.mtrand.RandomState(seed).get_state()
np.random.set_state(rs)
sigma = 0.1
ydata = ysim['A_obs'] + (np.random.randn(len(ysim['A_obs'])) * sigma)

# Get an instance of the ODE solver
solver = Solver(model, tspan)

# Variances of the datapoints for the likelihood calculation
variances = np.array([0.1] * len(tspan)) ** 2

# The likelihood function
def likelihood(x):
    """The argument x is a vector of parameter values at the current step."""
    # x[0] == A_0
    # x[1] == k
    # Transform the parameters back to a linear scale
    lin_x = 10 ** x
    solver.run([lin_x[0], lin_x[1], 1.0])
    # Calculate the likelihood
    return -np.sum(((ydata - solver.yobs['A_obs'])**2) / (2 * variances))

# Define the parameters of the walk
nwalkers = 10
ndim = 2
# Define the set of initial guesses
x0 = np.random.uniform(size=(nwalkers, ndim))
# Log-transform the guesses
x0 = np.log10(x0)

sampler = emcee.EnsembleSampler(nwalkers, ndim, likelihood)
sampler.random_state = rs

print "Burn in sampling..."
pos, prob, state = sampler.run_mcmc(x0, 10)
sampler.reset()

print "Main sampling..."
sampler.run_mcmc(x0, 100)

print "Done sampling."

# Shorthand variable
chain = 10 ** sampler.flatchain

# Plot histograms of the parameters
plt.ion()
plt.figure()
plt.subplot(1, 2, 1)
plt.hist(chain[:, 0], 100, color='k', histtype='step')
plt.xlabel('A_0')
plt.subplot(1, 2, 2)
plt.hist(chain[:, 1], 100, color='k', histtype='step')
plt.xlabel('k')

print("Mean acceptance fraction: {0:.3f}".format(
                    np.mean(sampler.acceptance_fraction)))

# Plot the original data (underlying and noisy)
# along with the sampled trajectories
plt.figure()
plt.plot(tspan, ysim, color='r')
plt.plot(tspan, ydata, color='g')

# Specify the number of sampled trajectories to plot
num_to_plot = 100
if num_to_plot > len(chain):
    num_to_plot = len(chain)

# Plot the sampled trajectories
for i in range(num_to_plot):
    # Get a random sample from the walk
    rand_ix = np.random.randint(len(chain))
    A0 = chain[rand_ix, 0]
    k = chain[rand_ix, 1]
    # Run the solver with the sampled parameters
    solver.run([A0, k, 1.0])
    # Plot alongside the data
    plt.plot(tspan, solver.yobs['A_obs'], alpha=0.05, color='b')

