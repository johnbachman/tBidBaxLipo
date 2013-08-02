from matplotlib import pyplot as plt
from tbidbaxlipo.simdata.sim_test.run_script import Job

# Create the job instance
j = Job()

# Run the deterministic simulation
(t, det_obs) = j.run_one_cpt()

# Plot deterministic results
plt.ion()
plt.figure()
plt.plot(t, det_obs['pores'], label='one_cpt')

# Run the stochastic simulation
xrecs = j.run_n_cpt(cleanup=True)
(means, stds) = j.calculate_mean_and_std(xrecs)

# Plot stochastic results
plt.errorbar(means['time'], means['pores'] / j.scaling_factor,
             yerr=stds['pores'] / j.scaling_factor, label='n_cpt')

# Label the plot
plt.xlabel('Time (secs)')
plt.ylabel('Total pores')
plt.title('Comparing one_cpt and n_cpt simulations')
plt.legend(loc='lower right')
