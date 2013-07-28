# Put each condition in directory name?
# Probably better to put conditions in a list or array
# Then send jobs with an index into the param array
# And then reassemble jo

# Number of simulations to run
xrecs = []
dr_all = []

Get directory of sim results, read *.out
    sim_result = kappa.parse_kasim_outfile(out_filename)
    xrecs.append(sim_result)

# Convert the multiple simulations in an array...
xall = array([x.tolist() for x in xrecs])

# ...and calculate the Mean and SD across the simulations
x_std = recarray(xrecs[0].shape, dtype=xrecs[0].dtype, buf=std(xall, 0))
x_avg = recarray(xrecs[0].shape, dtype=xrecs[0].dtype, buf=mean(xall, 0))


# Iterate
    # Call bsub on each one

# Parse data
# Take job name
    # Read each file
    # Pool data into a large numpy array


# Run over a large number of conditions
# Specify conditions
    # Iterate over conditions
    # Spawn job

