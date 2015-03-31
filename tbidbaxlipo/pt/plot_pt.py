import sys
import pickle
import numpy as np
from matplotlib import pyplot as plt

if __name__ == '__main__':
    # Check arguments
    if len(sys.argv) < 2:
        print "Usage: %s mcmc_filename"
        sys.exit(1)
    # Open the MCMC file
    mcmc_filename = sys.argv[1]
    (gf, samp) = pickle.load(open(mcmc_filename))
    # The chain
    c = samp.chain
    ntemps = samp.chain.shape[0]
    # Lists to store the mean and SD of the log likelihood at each beta
    lkl_m = []
    lkl_sd = []
    # Iterate over the temperatures
    for t_ix in range(ntemps):
        lkl_m.append(np.mean(samp.lnlikelihood[t_ix]))
        lkl_sd.append(np.std(samp.lnlikelihood[t_ix]))
    # Convert lists to arrays
    lkl_m = np.array(lkl_m)
    lkl_sd = np.array(lkl_sd)
    # Get the log evidence from the thermodynamic integral
    (evi_m, evi_err) = samp.thermodynamic_integration_log_evidence(fburnin=0)
    plt.ion()
    # Plot the mean log likelihood vs. beta
    plt.figure()
    plt.errorbar(np.log10(samp.betas), lkl_m, yerr=lkl_sd, marker='o')
    plt.xlabel('log10(beta)')
    plt.ylabel('E[ln(lkl)]')
    plt.title('TI curve: log(P(D|M)) = %f +/- %f' % (evi_m, evi_err))
    # Plot the log10 of the mean log likelihood vs. beta (only useful
    # for getting a sense of the trend
    plt.figure()
    plt.plot(np.log10(samp.betas), np.log10(-lkl_m), marker='o')
    plt.xlabel('log10(beta)')
    plt.ylabel('log10(E[ln(lkl)])')
    plt.title('TI temperature curve, log y-axis')
