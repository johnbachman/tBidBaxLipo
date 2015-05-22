from emcee import PTSampler
import sys
import os
import socket
from emcee.mpi_pool import MPIPool
import numpy as np
import timeit
import time

# Get the memory usage as a string, works on Linux machines (but not MacOS)
def get_memory(pid):
    with open('/proc/%s/status' % pid) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('VmSize'):
                vmsize = line.strip()
    return vmsize

# Add some detail to logging messages
# (Helps see whether all processes are running on same host or not)
def log(msg, procinfo):
    (hostname, pid, rank) = procinfo
    print("HOST %s: PID %s: RANK: %d: %s" % (hostname, pid, rank, msg))

# Dummy likelihood and prior functions
def lkl(p):
    time.sleep(0.001)
    return 0

def pr(p):
    return 0

# Define the size of the sampling job, and the initial position
ntemps = 50
nwalkers = 400
ndim = 10
p0 = np.zeros((ntemps, nwalkers, ndim))

# Create the pool, use the loadbalancing mode
pool = MPIPool(loadbalance=True)

# Get the identity of this process
hostname = socket.gethostname()
pid = os.getpid()
rank = pool.rank
procinfo = (hostname, pid, rank)

poolsize = pool.size + 1 # Emcee's MPIPool subtracts one, so we add it back
print("MPI pool size: %d" % poolsize)

if not pool.is_master():
    log("Waiting", procinfo)
    pool.wait()
    sys.exit(0)

# Run the expt three times (helps to see if memory grows with repeated runs)
times = []
for i in range(3):
    log("Waiting to start", procinfo)
    time.sleep(5)

    # Create the sampler
    log("Creating the sampler", procinfo)
    sampler = PTSampler(ntemps, nwalkers, ndim, lkl, pr, pool=pool)

    # Run the sampler
    cur_start_position = p0
    num_steps = 50
    nstep = 0
    log("Starting timer", procinfo)
    start_time = timeit.default_timer()
    for p, lnprob, lnlike in sampler.sample(cur_start_position,
                                        iterations=num_steps, storechain=True):
        print("Step %d" % nstep)
        nstep += 1

    # Print the time elapsed
    elapsed = timeit.default_timer() - start_time
    times.append(elapsed)
    log("Time elapsed: %s" % elapsed, procinfo)

    # Print the memory
    log(get_memory(pid), procinfo)

# Write the results to a file
with open('%d.txt' % poolsize, 'w') as f:
    for t in times:
        f.write('%s\n' % t)

# Close the pool so that the process doesn't hang
pool.close()
sys.exit(0)

