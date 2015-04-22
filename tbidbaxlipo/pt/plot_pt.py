import sys
import pickle
import numpy as np
from matplotlib import pyplot as plt
from scipy.special import erfc
from scipy.interpolate import interp1d

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
    log_betas = np.log10(samp.betas)
    p_acc = samp.tswap_acceptance_fraction

    # Calculate average final posterior across chains
    # Mean probability over final points of all chains
    #final_probs = np.mean(samp.lnprobability[:, :, -1], axis=1)
    #final_prob_diffs = np.diff(final_probs)
    # Calculate standard deviation of the log probability
    #final_prob_sd = np.std(samp.lnprobability[:, :, -1], axis=1)
    #plt.figure()
    #plt.errorbar(log_betas, final_probs, yerr=final_prob_sd)
    E = -(samp.lnprobability[:, :, -1] / samp.betas[:, None])
    E_m = np.mean(E, axis=1)
    E_sd = np.std(E, axis=1)

    # Create interpolation functions for E(T) and sigma(T)
    E_B = interp1d(log_betas[::-1], E_m[::-1])
    sigma_B = interp1d(log_betas[::-1], E_sd[::-1])

    beta_grid = np.linspace(0, -6, 6e3)
    optimal_betas = []
    # Target swap acceptance rate of 0.3 implies dE/sigma of about 1.5
    # Start with energy at log beta of 0 (beta of 1)
    target_dE_sigma = 2.0
    # Initialize the beta, energy, and sigma values for the first point
    # (lowest temperature)
    cur_B = beta_grid[0]
    optimal_betas.append(cur_B)
    cur_E = E_B(cur_B) # Energy at the current beta
    cur_sigma = sigma_B(cur_B) # SD of energy at the current beta
    # For every point that we sample in the beta grid (skipping the first
    # point, since it's our reference)...
    diffs = []
    for cur_B in beta_grid[1:]:
        # ...calculate the new energy...
        new_E = E_B(cur_B)
        # ...and the new energy SD at this beta value.
        new_sigma = sigma_B(cur_B)
        # Calculate the mean sigma for the previous and current betas
        sigma_m = (new_sigma + cur_sigma) / 2.0 # FIXME
        # Now use these values to get the dE / sigma_m value at the proposed
        # beta point:
        dE_sigma = (new_E - cur_E) / sigma_m
        # Presumably the dE_sigma value will be small if the temperatures are
        # very close together, so we just check if the dE_sigma value is
        # greater than our target; if so, we update the beta and energy values.
        if dE_sigma > target_dE_sigma:
            cur_E = new_E
            cur_sigma = new_sigma
            optimal_betas.append(cur_B)
            diffs.append(dE_sigma)
    # Finally, add the last beta point in the grid to make sure we get full
    # coverage of the temperature range
    optimal_betas.append(beta_grid[-1])
    optimal_betas = np.array(optimal_betas)
    #with open('optimal_betas.pck', 'w') as f:
        #pickle.dump(10 ** optimal_betas, f)
    diffs = np.array(diffs)

    # Plot the log10 of the mean log likelihood vs. beta (only useful
    # for getting a sense of the trend)
    plt.figure()
    logloglkl = np.log10(-lkl_m)
    lkl_B = interp1d(log_betas[::-1], logloglkl[::-1])
    plt.plot(log_betas, logloglkl, marker='.')
    plt.plot(optimal_betas, lkl_B(optimal_betas), color='red', marker='o')
    plt.xlabel('log10(beta)')
    plt.ylabel('log10(E[ln(lkl)])')
    plt.title('TI temperature curve, log y-axis')

    # Calculate the area of overlap
    num_diffs = len(E) - 1
    overlap = np.zeros(num_diffs)
    sigma_m = np.zeros(num_diffs)
    dE = np.zeros(num_diffs)
    for i in xrange(num_diffs):
        dE[i] = E_m[i + 1] - E_m[i]
        sigma_m[i] = (E_sd[i] + E_sd[i+1]) / 2.0
        overlap[i] = erfc(dE[i] / (2.0 * np.sqrt(2) * sigma_m[i]))

    plt.figure()
    plt.plot(log_betas, E_m, marker='o')
    plt.title('E(beta)')

    plt.figure()
    plt.plot(log_betas, E_sd, marker='o')
    plt.title('sigma(beta)')

    plt.figure()
    plt.plot(overlap, p_acc[:-1], linestyle='', marker='o')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xlabel('overlap')
    plt.ylabel('p_acc')
    plt.title('acceptance vs. overlap')

    plt.figure()
    plt.plot(dE / sigma_m, p_acc[:-1], marker='o')
    plt.ylim([0, 1])
    plt.xlabel('dE/sigma')
    plt.ylabel('p_acc')
    plt.title('acceptance vs. dE/sigma')

    # Plot the mean log likelihood vs. beta
    plt.figure()
    plt.errorbar(log_betas, lkl_m, yerr=lkl_sd, marker='o')
    plt.xlabel('log10(beta)')
    plt.ylabel('E[ln(lkl)]')
    plt.title('TI curve: log(P(D|M)) = %f +/- %f' % (evi_m, evi_err))

    # Plot the log10 of the mean log likelihood vs. beta (only useful
    # for getting a sense of the trend
    plt.figure()
    plt.plot(log_betas, np.log10(-lkl_m), marker='o')
    plt.xlabel('log10(beta)')
    plt.ylabel('log10(E[ln(lkl)])')
    plt.title('TI temperature curve, log y-axis')

    # Calculate all differences between betas
    beta_diffs = np.diff(log_betas)
    lkl_diffs = np.diff(lkl_m * samp.betas)
    slopes = lkl_diffs / beta_diffs
    plt.figure()
    plt.plot(slopes, samp.tswap_acceptance_fraction[:-1], linestyle='',
             marker='o')

    # Calculate all differences between log(log(lkl))
    log_lkl_diffs = np.diff(np.log10(-lkl_m))
    slopes = log_lkl_diffs / beta_diffs
    plt.figure()
    plt.plot(slopes, samp.tswap_acceptance_fraction[:-1], marker='o')
    plt.title('loglkl/beta diffs')

