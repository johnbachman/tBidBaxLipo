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
        plt.ion()
        file_basename = None

    # Plot all (raw data)
    if plot_type == 'raw':
        nba.plot_all(df, nbd_residues, file_basename=file_basename)
    # Plot endpdoints
    elif plot_type == 'endpoints':
        nba.plot_nbd_endpoints(df, nbd_residues, last_n_pts=3,
                           file_basename=file_basename)
        nba.plot_release_endpoints(df, nbd_residues, normalized_to_wt=False,
                           last_n_pts=3, file_basename=file_basename)
    # Unknown plot type
    else:
        print "Unknown plot type: %s" % plot_type
        sys.exit(1)
