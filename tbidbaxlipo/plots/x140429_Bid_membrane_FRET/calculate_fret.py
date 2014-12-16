from tbidbaxlipo.data.parse_140429_Bid_membrane_FRET import \
        time, fda, fa, bg, fd, timecourse_wells
import itertools
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.util.plate_assay import *
from tbidbaxlipo.util import fitting
from tbidbaxlipo.util.error_propagation import calc_ratio_mean_sd

### PLOTTING FUNCTIONS ### 

def plot_fd_final_values(wells, fd):
    """Plot the endpoint fluorescence of the FD (donor-only) condition, to show
    that it fluctuates non-monotonically as a function of the position of the
    well in the dilution series (i.e., as a function of WT Bid
    concentration).
    """
    (fd_avgs, fd_stds) = averages(wells, fd)
    bid_concs =  []
    final_val_avgs = []
    for conc_name in fd_avgs:
        bid_concs.append(float(conc_name.split(' ')[1]))
        final_val_avg = np.mean(fd_avgs[conc_name][VALUE][-10:])
        final_val_avgs.append(final_val_avg)
    plt.figure('FD values, final')
    plt.plot(bid_concs, final_val_avgs, marker='o')
    plt.xlabel('[WT Bid] (nM)')
    plt.ylabel('RFU')
    plt.xlim([-10, 1010])

def get_average_fa(wells, fa, plot=True):
    """Plot all FA (acceptor only, no donor) trajectories and then return the
    average across all of them. This average is used to do background subtraction
    for the FDA (donor + acceptor) trajectories.

    The fluorescence of the FA condition is not sensitive to the amount of
    unlabeled cBid present, so we can calculate our background vector by
    averaging across all of the conditions.

    Parameters
    ----------
    wells : collections.OrderedDict
        Dict containing the timecourses in each well.
    fa : collections.OrderedDict
        Dict mapping the concentration conditions to sets of wells for
        the FA condition.
    plot : boolean
        Whether to plot the timecourses.

    Return
    ------
    numpy.array containing the average FA timecourse.
    """

    # Calculate average background
    fa_well_names = itertools.chain(*fa.values())
    fa_wells = extract(fa_well_names, timecourse_wells)
    fa_avg = np.zeros((len(fa_wells.keys()), len(time)))
    for i, well_name in enumerate(fa_wells.keys()):
        fa_avg[i] = wells[well_name][VALUE]
    fa_avg = np.mean(fa_avg, axis=0)
    # Plot
    if plot:
        plt.figure('FA')
        plot_all(fa_wells)
        plt.plot(time, fa_avg, linewidth=3, color='k')
    # Return
    return fa_avg

def get_average_bg(wells, bg, plot=True):
    """Plot all BG (no donor, no acceptor) trajectories and then return the
    average across all of them. This average is used to do background subtraction
    for the FD (donor only, no acceptor) trajectories.

    The fluorescence of the BG condition is not sensitive to the amount of
    unlabeled cBid present, so we can calculate our background vector by
    averaging across all of the conditions.

    Parameters
    ----------
    wells : collections.OrderedDict
        Dict containing the timecourses in each well.
    bg : collections.OrderedDict
        Dict mapping the concentration conditions to sets of wells for
        the FA condition.
    plot : boolean
        Whether to plot the timecourses.

    Return
    ------
    numpy.array containing the average BG timecourse.
    """
    # Calculate average background
    bg_well_names = itertools.chain(*bg.values())
    bg_wells = extract(bg_well_names, timecourse_wells)
    bg_avg = np.zeros((len(bg_wells.keys()), len(time)))
    for i, well_name in enumerate(bg_wells.keys()):
        bg_avg[i] = wells[well_name][VALUE]
    bg_avg = np.mean(bg_avg, axis=0)
    # Plot
    if plot:
        plt.figure('BG')
        plot_all(bg_wells)
        plt.plot(time, bg_avg, linewidth=3, color='k')
    # Return
    return bg_avg

def calculate_fret_from_endpoints(wells, fda, fd, fa, bg, num_pts=20,
                                  plot=True):
    """Calculate FRET using the endpoints of the FDA and FD conditions.

    For the donor + acceptor, subtracts the acceptor-only background average;
    for the donor-only, subtracts the no donor/no acceptor background average.
    Then takes the average of the last num_pts of the background-subtracted traces
    to get the fluorescence value for that condition/replicate.

    To get an estimate of error across replicates, the values for the FDA and
    FA wells are then averaged across replicates, yielding a mean fluorescence
    for each with an associated standard error. The FRET is then calculated by
    taking the ratio, and the standard error of the FRET value is calculated by
    sampling the standard deviation of the ratio of two normal distributions
    with means and standard deviations matching those of the mean and standard
    errors of the FDA and FD conditions.

    Note that this approach accounts for error across replicates, but not for
    error in determining the mean fluorescence of each individual replicate.
    However, it appears that the variability between replicates is
    substantially greater than the variability within a single fluorescence
    timecourse, and the standard error of the mean could be further reduced by
    averaging across more points.

    Parameters
    ----------
    wells : collections.OrderedDict
        Dict containing the timecourses for each well in the plate.
    fda, fd, fa, bg : collections.OrderedDict
        Dicts mapping the concentration conditions to the well names of the
        replicates.
    num_pts : int
        The number of final points to average across.
    plot : boolean
        Whether to plot the FRET results (defaults to True)

    Returns
    -------
    tuple : (bid_concs, fret_means, fret_ses)
        Tuple containing the concentrations of the unlabeled Bid used in each
        conditions; the mean FRET values calculated for each unabeled Bid
        concentration; and the standard error of the mean calculated by
        sampling.
    """

    # Get the background averages
    bg_avg = get_average_bg(wells, bg, plot=False)
    fa_avg = get_average_fa(wells, fa, plot=False)
    # Initialize a few variables
    num_concs = len(fda.keys())
    num_reps = 3
    bid568_conc = 10. # Concentration of Bid-568, in nM
    fda_endpoints = np.zeros((num_concs, num_reps))
    fd_endpoints = np.zeros((num_concs, num_reps))
    bid_concs = np.zeros(num_concs)
    # Iterate over the conditions
    for conc_ix, conc_name in enumerate(fda.keys()):
        bid_concs[conc_ix] = float(conc_name.split(' ')[1])
        # FDA
        for rep_ix, well_name in enumerate(fda[conc_name]):
            fda_value = np.array(wells[well_name][VALUE])
            # Subtract the FA background:
            fda_value = fda_value - fa_avg
            # Calculate endpoint value over last num_pts
            endpt = np.mean(fda_value[-num_pts:])
            fda_endpoints[conc_ix, rep_ix] = endpt
        # FD
        for rep_ix, well_name in enumerate(fd[conc_name]):
            fd_value = np.array(wells[well_name][VALUE])
            # Subtract the BG background:
            fd_value = fd_value - bg_avg
            # Calculate endpoint value over last num_pts
            endpt = np.mean(fd_value[-num_pts:])
            fd_endpoints[conc_ix, rep_ix] = endpt

    # Now, calculate the means of the endpts across replicates
    fda_means = np.mean(fda_endpoints, axis=1)
    fd_means = np.mean(fd_endpoints, axis=1)
    fda_ses = np.std(fda_endpoints, axis=1, ddof=1) / np.sqrt(float(num_reps))
    fd_ses = np.std(fd_endpoints, axis=1, ddof=1) / np.sqrt(float(num_reps))
    # Calculate the mean and SEM of the FDA/FD ratio distributions
    fret_tuples = [calc_ratio_mean_sd(fda_means[i], fda_ses[i],
                                      fd_means[i], fd_ses[i])
                   for i in range(len(fda_means))]
    # Convert the list of tuples into two separate lists
    (fret_means, fret_ses) = zip(*fret_tuples)
    # Now convert to numpy arrays. We take 1 - the FRET values so that we get
    # FRET in terms of FRET efficiency instead of % quenching:
    fret_means = 1 - np.array(fret_means)
    fret_ses = np.array(fret_ses)

    if plot:
        plt.figure('FDA and FD, mean/SEM over reps')
        plt.errorbar(np.log10(bid_concs), fda_means, yerr=fda_ses)
        plt.errorbar(np.log10(bid_concs), fd_means, yerr=fd_ses)
        plt.xlabel('log10([WT Bid]) (nM)')
        plt.ylabel('RFU')
        plt.xlim([-0.5, 3.1])
        plt.title('FDA and FD, mean/SEM over reps')

        plt.figure('FRET')
        plt.errorbar(np.log10(bid_concs), fret_means, yerr=fret_ses)
        plt.xlabel('log10([WT Bid]) (nM)')
        plt.ylabel('FRET (%)')
        plt.xlim([-0.5, 3.1])
        plt.title('FRET')

    return (bid_concs, fret_means, fret_ses)

def fit_fda_and_fd_curves(wells, fda, fd, fa, bg):
    # Get the background averages
    bg_avg = get_average_bg(wells, bg, plot=False)
    fa_avg = get_average_fa(wells, fa, plot=False)
    # Now, plot the curves, with fits
    for conc_name in fda.keys():
        plt.figure()
        # FDA
        for well_name in fda[conc_name]:
            fda_value = np.array(wells[well_name][VALUE])
            # Subtract the FA background:
            fda_value = fda_value - fa_avg
            # Do the fit
            f0 = fitting.Parameter(3500.)
            fmax = fitting.Parameter(2000.)
            k = fitting.Parameter(1e-3)
            def fit_func(t):
                return fmax()*np.exp(-k()*t) + f0()
            res = fitting.fit(fit_func, [f0, fmax, k], fda_value, time)
            # Plot data and fit
            plt.plot(time, fda_value)
            plt.plot(time, fit_func(time), color='k', linewidth=2)
        # FD
        for well_name in fd[conc_name]:
            fd_value = np.array(wells[well_name][VALUE])
            # Subtract the BG background:
            fd_value = fd_value - bg_avg
            # Do the fit
            f0 = fitting.Parameter(3500.)
            fmax = fitting.Parameter(2000.)
            k = fitting.Parameter(1e-3)
            def fit_func(t):
                return fmax()*np.exp(-k()*t) + f0()
            res = fitting.fit(fit_func, [f0, fmax, k], fd_value, time)
            # Plot data and fit
            plt.plot(time, fd_value)
            plt.plot(time, fit_func(time), color='k', linewidth=2)

def get_fret_from_endpoints(num_pts=20):
    """Wrapping the call to calculate_fret_from_endpoints allows us to get
    the FRET data from this module without explicitly providing all of the
    necessary arguments to the function."""
    return calculate_fret_from_endpoints(timecourse_wells, fda, fd, fa, bg,
                                         num_pts=num_pts, plot=False)

if __name__ == '__main__':
    plt.ion()

    plot_fd_final_values(timecourse_wells, fd)
    bg_avg = get_average_bg(wells, bg, plot=False)
    fa_avg = get_average_fa(wells, fa, plot=False)

    calculate_fret_from_endpoints(timecourse_wells, fda, fd, fa, bg, num_pts=20,
                                  plot=True)

