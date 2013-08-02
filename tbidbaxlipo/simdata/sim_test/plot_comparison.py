from matplotlib import pyplot as plt
from tbidbaxlipo.simdata.sim_test import means, stds
from tbidbaxlipo.simdata.sim_test.run_script import Job

# Create the job instance
j = Job()

# Run the deterministic simulation
(t, det_obs) = j.run_one_cpt()

# Plot deterministic results
plt.ion()
plt.figure()
plt.plot(t, det_obs['pores'])

# Plot stochastic results
plt.errorbar(means['time'], means['pores'] / j.scaling_factor,
             yerr=stds['pores'] / j.scaling_factor)

