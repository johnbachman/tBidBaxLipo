import sys
from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
import tbidbaxlipo.plots.nbd_bax_analysis as nba
from matplotlib import pyplot as plt

if __name__ == '__main__':
    # Usage
    usage = "Usage: python plots.py [plot_type] [(optional) file_basename]"
    if len(sys.argv) < 2:
        print usage
        sys.exit(1)

    # Which type of plot to make?
    plot_type = sys.argv[1]

    # Set the file basename; if not specified, assume we're plotting
    # interactively
    if len(sys.argv) > 2:
        file_basename = sys.argv[2]
    else:
        file_basename = None

    # Plot all (raw data)
    if plot_type == 'raw':
        nba.plot_all(df, nbd_residues, file_basename=file_basename)
    # Plot endpoints
    elif plot_type == 'nbd_endpoint_norm':
        nba.plot_nbd_endpoints(df, nbd_residues, last_n_pts=3,
                           file_basename=file_basename, normalize_nbd=True)
    elif plot_type == 'nbd_endpoint_no_norm':
        nba.plot_nbd_endpoints(df, nbd_residues, last_n_pts=3,
                           file_basename=file_basename, normalize_nbd=False)
    elif plot_type == 'release_endpoint_norm':
        nba.plot_release_endpoints(df, nbd_residues, normalized_to_wt=True,
                           last_n_pts=3, file_basename=file_basename)
    elif plot_type == 'release_endpoint_no_norm':
        nba.plot_release_endpoints(df, nbd_residues, normalized_to_wt=False,
                           last_n_pts=3, file_basename=file_basename)
    elif plot_type == 'initial_rate_samples':
        nba.plot_initial_rate_samples(df, nbd_residues, timepoint_ix=4,
                                      file_basename=file_basename)
    elif plot_type == 'calc_release_peaks':
        nba.calc_release_peaks(df, nbd_residues,
                               csv_filename='data1_release_peak_times.csv')
    elif plot_type == 'example_derivatives':
        nba.plot_example_derivatives(df, 'Bid', '15', 1,
                                     plot_filename='data1_derivatives_Bid_15_r1')
        nba.plot_example_derivatives(df, 'Bid', '54', 1,
                                     plot_filename='data1_derivatives_Bid_54_r1',
                                     plot_tb_peak=True)
    elif plot_type == 'evidence':
        from tbidbaxlipo.plots.bid_bim_nbd_release.plot_bf_values import \
             plot_2confs, plot_3confs
        plot_2confs('data1_evidence_2confs')
        plot_3confs('data1_evidence_3confs')
    # Unknown plot type
    else:
        print "Unknown plot type: %s" % plot_type
        sys.exit(1)
