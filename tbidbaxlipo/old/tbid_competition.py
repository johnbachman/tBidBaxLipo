from pysb import *
import numpy as np
from matplotlib import pyplot as plt
from pysb.integrate import Solver

# Use interactive plotting mode
plt.ion()

# Create the pysb model
Model()

# Two entities: Bid (labeled and unlabeled) and binding sites for Bid
Monomer('Bid', ['a568', 'bs', 'nonsp'], {'a568': ['y', 'n'], 'nonsp':['y','n']})
Monomer('Bid_site', ['bs'])
Monomer('Liposome', [])


# Parameters defining the initial concentrations, in nM
Parameter('Bid_568_0', 10)
Parameter('Bid_unlab_0', 0)
Parameter('Liposomes_0', 1)
Parameter('Sites_per_lipo', 40)
Expression('Bid_site_0', Liposomes_0 * Sites_per_lipo)

# Define the initial conditions using the above parameters
Initial(Bid(a568='y', bs=None, nonsp='n'), Bid_568_0)   # labeled Bid
Initial(Bid(a568='n', bs=None, nonsp='n'), Bid_unlab_0) # unlabeled Bid
Initial(Liposome(), Liposomes_0)
Initial(Bid_site(bs=None), Bid_site_0)       # Bid binding sites

# Another possibility: nonspecific biding to lipos
Parameter('Bid_binding_nonsp_kf', 1e-2) # Diffusion limited on-rate
Parameter('Bid_binding_nonsp_kr', 1e-2) # Avg lifetime of 100 sec for complex
# Our only reaction rule: Bid (whether labeled or unlabeled) binds to sites
Rule('Bid_binds_nonsp_kf',
     Bid(bs=None, nonsp='n') + Liposome() >>
     Bid(bs=None, nonsp='y') + Liposome(),
     Bid_binding_nonsp_kf)
Rule('Bid_binds_nonsp_kr',
     Bid(bs=None, nonsp='y') >> Bid(bs=None, nonsp='n'),
     Bid_binding_nonsp_kr)


Parameter('Bid_binding_kf', 1e-3) # Diffusion limited on-rate
Parameter('Bid_binding_kr', 1e-4) # Avg lifetime of 100 sec for complex
# Our only reaction rule: Bid (whether labeled or unlabeled) binds to sites
Rule('Bid_binds_sites',
     Bid(bs=None, nonsp='n') + Bid_site(bs=None) <>
     Bid(bs=1, nonsp='n') % Bid_site(bs=1),
     Bid_binding_kf, Bid_binding_kr)

# We only care about the fraction of the labeled cBid bound, so it's our
# only observable
Observable('Bid_568_bound', Bid(a568='y', bs=1, nonsp='n') % Bid_site(bs=1))
Observable('Bid_568_nonsp_bound', Bid(a568='y', bs=None, nonsp='y'))
Expression('Bid_FRET', Bid_568_bound + Bid_568_nonsp_bound)

# Define the concentrations of the competitor
unlab_Bid_dilution_series = np.linspace(0, 10000, 100)
#unlab_Bid_dilution_series[0] = 10000.
#unlab_Bid_dilution_series[11] = 0.
#for i in range(1, 11):
#    unlab_Bid_dilution_series[i] = unlab_Bid_dilution_series[i - 1] / 2.

# For competitor concentration, run the simulation and plot it
max_timepoint = 7200. # Duration of simulation, in seconds
t = np.linspace(0, 7200, 1000) # The time vector for the simulation
s = Solver(model, t)

# Plot the simulated binding timecourses
plt.figure()
endpoints = np.zeros(len(unlab_Bid_dilution_series))
for i, unlab_cBid_conc in enumerate(unlab_Bid_dilution_series):
    # Set the unlabeled cBid concentration in the model
    Bid_unlab_0.value = unlab_cBid_conc
    s.run()
    Bid_fret_timecourse = s.yexpr['Bid_FRET'] / Bid_568_0.value
    plt.plot(t, Bid_fret_timecourse)
    endpoints[i] = Bid_fret_timecourse[-1]
plt.title('Bid competition timecourses')

# Plot on a linear concentration scale
plt.figure()
plt.plot(unlab_Bid_dilution_series, endpoints, marker='o')
plt.title('Bid competition, linear scale')

# Now plot on a log scale (skip the 0 point because you can't take the log of 0)
plt.figure()
plt.plot(np.log10(unlab_Bid_dilution_series[:-1]), endpoints[:-1], marker='o')
plt.title('Bid competition, log scale')

# Do a simulated liposome titration experiment, with no competitor
num_lipo_concs = 50
#lipo_concs = np.zeros(num_lipo_concs)
#lipo_concs[0] = 30.
#lipo_concs[num_lipo_concs - 1] = 0.
#for i in range(1, num_lipo_concs - 1):
#    lipo_concs[i] = lipo_concs[i - 1] / 2.
lipo_concs = np.linspace(0, 30, num_lipo_concs)

# Plot the simulated binding timecourses
plt.figure()
Bid_endpoints = np.zeros(len(lipo_concs))
for i, lipo_conc in enumerate(lipo_concs):
    # Set the unlabeled cBid concentration in the model
    Bid_unlab_0.value = 0
    Liposomes_0.value = lipo_conc
    s.run()
    Bid_fret_timecourse = s.yexpr['Bid_FRET']
    Bid_fret_timecourse /= Bid_568_0.value
    plt.plot(t, Bid_fret_timecourse)
    Bid_endpoints[i] = Bid_fret_timecourse[-1]
plt.title('Bid liposome titration timecourses')

plt.figure()
plt.plot(lipo_concs, Bid_endpoints, marker='o')


