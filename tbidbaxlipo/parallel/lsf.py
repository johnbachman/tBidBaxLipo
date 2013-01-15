"""Organize the submission of multiple calibration jobs to LSF, and recombine
the resulting data back.

The MCMC object should contain everything necessary to re-run the fitting
procedure if necessary (not just the steps/positions, but also the model,
initial parameter values, perhaps even the data?).
"""

def submit(queue='sysbio_unlimited'):
   pass 


def multi_start_fits(n_fits=10):
    for fit_num in range(0, n_fits):
        # Tell do_fit to pick a random starting point
        pass

    # For multiple j

