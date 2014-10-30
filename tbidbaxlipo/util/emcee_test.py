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
#np.random.rs = np.random.mtrand.RandomState(seed).get_state()
#np.random.set_state(rs)
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

def run_sampler(num_steps, random_state, pos=None):
    # Define the parameters of the walk
    nwalkers = 10
    ndim = 2

    if pos is None:
        # Initialize the random state here, with a seed? Or pass in the random
        # state
        # Define the set of initial guesses
        x0 = np.random.uniform(size=(nwalkers, ndim))
        # Log-transform the guesses
        x0 = np.log10(x0)
    else:
        x0 = pos

    sampler = emcee.EnsembleSampler(nwalkers, ndim, likelihood)
    sampler.random_state = random_state

    #print "Burn in sampling..."
    #pos, prob, state = sampler.run_mcmc(x0, 10)
    #sampler.reset()

    print "Main sampling..."
    sampler.run_mcmc(x0, num_steps)

    print "Done sampling."
    return sampler


if __name__ == '__main__':
    # This is what happens when you seed with 1 and then run for 100 steps.
    seed = 2
    np.random.seed(seed)
    rs = np.random.get_state()

    num_steps = 100
    sampler1 = run_sampler(num_steps, rs)
    # Final position
    x49_1 = sampler1.chain[:, 49]
    x99_1 = sampler1.chain[:, 99]

    # Now let's see what happens when you seed with 1, run for 50 steps, then
    # reseed to continue:
    np.random.seed(seed)
    rs = np.random.get_state()
    num_steps = 50
    sampler2 = run_sampler(num_steps, rs)
    x49_2 = sampler2.chain[:, 49]

    #print x49_1
    #print x49_2
    #print x49_1 - x49_2

    # Now continue sampling the second half
    rs = sampler2.random_state
    sampler3 = run_sampler(num_steps, rs, pos=x49_2)
    x99_2 = sampler3.chain[:, 49]

    #print x99_1
    #print x99_2
    #print x99_1 - x99_2

    # Compare different trajectories
    plt.ion()
    plt.plot(np.arange(100), sampler1.chain[:,:,0].T, color='r')
    plt.plot(np.arange(50), sampler2.chain[:,:,0].T, color='g',
            linestyle='--')
    plt.plot(np.arange(50, 100), sampler3.chain[:,:,0].T, color='b',
            linestyle='--')
    plt.title('Comparison of continuous vs. interrupted chains')

    # Shorthand variable
    chain = 10 ** sampler1.flatchain

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
                        np.mean(sampler1.acceptance_fraction)))

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

