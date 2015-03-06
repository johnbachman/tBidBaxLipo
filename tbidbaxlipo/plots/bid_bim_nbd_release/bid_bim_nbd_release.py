from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import fitting
from tbidbaxlipo.models.nbd.multiconf import Builder
from pysb.integrate import Solver
import itertools
import titration_fits as tf
import scipy.stats
import scipy.signal
from tbidbaxlipo.util.error_propagation import calc_ratio_sd, calc_ratio_mean_sd
from tbidbaxlipo.util import moving_average
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis, \
                             fontsize

line_colors = {'Bid': 'r', 'Bim': 'b'}
line_styles = {1:':', 2:'-', 3:'--'}
rep_colors = {1:'r', 2:'g', 3:'b'}

def plot_derivatives(activator):
    nbd_sites = ['WT', '3', '5', '15', '36', '40', '47', '54', '62', '68',
                 '79', '120', '122', '126', '138', '151', '175', '179',
                 '184', '188']
    replicates = range(1, 4)
    num_pts = 4
    window = 5 # For moving average
    for nbd_index, nbd_site in enumerate(nbd_sites):
        rn_ratios = []
        r_maxs = []
        n_maxs = []
        #rn_errs = []
        for rep_index in replicates:
            rt = df[(activator, 'Release', nbd_site, rep_index, 'TIME')].values
            ry = df[(activator, 'Release', nbd_site, rep_index, 'VALUE')].values
            # Take the moving average of the timecourse
            r_avg = moving_average(ry, n=window)
            # Take the derivative
            r_diff = np.diff(r_avg)
            # Calculate the maximum rate of change
            r_max = np.max(np.abs(r_diff))
            r_maxs.append(r_max)

            # Calculate max NBD slope, but not for WT
            if nbd_site == 'WT':
                n_maxs.append(0)
                rn_ratios.append(0)
            else:
                nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
                ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values
                n_avg = moving_average(ny, n=window)
                n_diff = np.diff(n_avg)
                n_max = np.max(np.abs(n_diff))
                n_maxs.append(n_max)

                rn_ratio = r_max / n_max
                rn_ratios.append(rn_ratio)

            plt.figure('%s, NBD-%s-Bax derivative' % (activator, nbd_site),
                       figsize=(12, 5))
            plt.subplot(1, 2, 1)
            plt.plot(rt[1+window-1:], r_diff, color=rep_colors[rep_index],
                     label='%s Rep %d' % (activator, rep_index))
            plt.ylabel('dRel/dt (% rel $sec^{-1}$)')
            plt.title('NBD-%s-Bax, Tb derivative' % nbd_site)
            plt.legend(loc='upper right')

            if nbd_site != 'WT':
                plt.subplot(1, 2, 2)
                plt.plot(nt[1+window-1:], n_diff, color=rep_colors[rep_index],
                         label='%s Rep %d' % (activator, rep_index))
                plt.xlabel('Time (sec)')
                plt.ylabel('dNBD/dt ($F/F_0\ sec^{-1}$)')
                plt.title('NBD-%s-Bax, NBD derivative' % nbd_site)
                plt.legend(loc='upper right')

                # Plot normalized slopes
                plt.figure('%s, NBD-%s-Bax normalized derivative' %
                           (activator, nbd_site))
                n_diff_norm = n_diff / np.max(np.abs(n_diff))
                r_diff_norm = r_diff / np.max(np.abs(r_diff))
                plt.plot(rt[1+window-1:], r_diff_norm,
                         color=rep_colors[rep_index],
                         linestyle=line_styles[2],
                         label='%s Rep %d' % (activator, rep_index))
                plt.plot(nt[1+window-1:], n_diff_norm,
                         color=rep_colors[rep_index], linestyle=line_styles[3])
                plt.xlabel('Time (sec)')
                plt.ylabel('% max rate')
                plt.title('%s, NBD-%s-Bax normalized derivative' %
                          (activator, nbd_site))
                plt.legend(loc='upper right')
        plt.figure('%s, NBD-%s-Bax derivative' % (activator, nbd_site))
        plt.tight_layout()

        plt.figure("Tb/NBD peak slope ratio")
        if activator == 'Bid':
            plt.bar(range(nbd_index*7, (nbd_index*7) + 3), rn_ratios,
                    width=1, color='r')
                    #yerr=rn_errs,
        else:
            plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), rn_ratios,
                    width=1, color='g')
                    #yerr=rn_errs,
        plt.figure("Tb peak slope")
        if activator == 'Bid':
            plt.bar(range(nbd_index*7, (nbd_index*7) + 3), r_maxs,
                    width=1, color='r')
        else:
            plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), r_maxs,
                    width=1, color='g')
        plt.figure("NBD peak slope")
        if activator == 'Bid':
            plt.bar(range(nbd_index*7, (nbd_index*7) + 3), n_maxs,
                    width=1, color='r')
        else:
            plt.bar(range(nbd_index*7+3, (nbd_index*7) + 6), n_maxs,
                    width=1, color='g')
    num_sites = len(nbd_sites)
    fig_names = ["Tb/NBD peak slope ratio", "Tb peak slope", "NBD peak slope"]
    for fig_name in fig_names:
        plt.figure(fig_name)
        ax = plt.gca()
        ax.set_xticks(np.arange(3, 3 + num_sites * 7, 7))
        ax.set_xticklabels(nbd_sites)

def plot_filtered_data(activator):
    nbd_sites = ['15', '54']
    #nbd_sites = ['WT', '3', '5', '15', '36', '40', '47', '54', '62', '68',
    #             '79', '120', '122', '126', '138', '151', '175', '179',
    #             '184', '188']
    replicates = range(1, 4)
    num_pts = 4
    window = 5 # For moving average
    # Create an order 3 lowpass butterworth filter.
    b, a = scipy.signal.butter(1, 0.2)

    for nbd_index, nbd_site in enumerate(nbd_sites):
        for rep_index in replicates:
            """
            rt = df[(activator, 'Release', nbd_site, rep_index, 'TIME')].values
            ry = df[(activator, 'Release', nbd_site, rep_index, 'VALUE')].values
            r_filt = scipy.signal.filtfilt(b, a, ry)
            plt.figure()
            plt.plot(rt, ry, linestyle='', marker='.')
            plt.plot(rt, r_filt)
            """
            nt = df[(activator, 'NBD', nbd_site, rep_index, 'TIME')].values
            ny = df[(activator, 'NBD', nbd_site, rep_index, 'VALUE')].values
            n_filt = scipy.signal.filtfilt(b, a, ny)
            plt.figure()
            plt.plot(nt, ny, linestyle='', marker='.')
            plt.plot(nt, n_filt)

if __name__ == '__main__':
    plt.ion()
    #plot_all(df, nbd_residues, file_basename='data1_raw')
    #plot_endpoints(df, nbd_residues, file_basename='data1')
    #plot_initial_rates('Bid')
    #plot_initial_rates('Bim')
    #plot_filtered_data('Bid')
    #fr = plot_3conf_fits('Bid')
    #plot_3conf_fits('Bid')
    #plot_peak_slopes('Bim')
    plot_derivatives('Bid')
