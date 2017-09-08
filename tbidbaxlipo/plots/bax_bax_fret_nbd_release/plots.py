import sys
from tbidbaxlipo.plots.bax_bax_fret_nbd_release.preprocess_data \
        import df, df_pre, nbd_residues, activators
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
    # Plot endpoints
    elif plot_type == 'nbd_endpoint_norm':
        nba.plot_nbd_endpoints(df, nbd_residues, last_n_pts=3,
                           file_basename=file_basename, normalize_nbd=True,
                           activators=activators)
    # Not done b/c data is already normalized in spreadsheet
    #elif plot_type == 'nbd_endpoint_no_norm':
    #    nba.plot_nbd_endpoints(df, nbd_residues, last_n_pts=3,
    #                       file_basename=file_basename, normalize_nbd=False,
    #                       activators=activators)
    # Can't do normalized release because no WT reference in the dataset
    #elif plot_type == 'release_endpoint_norm':
    #    nba.plot_release_endpoints(df, nbd_residues, normalized_to_wt=True,
    #                       last_n_pts=3, file_basename=file_basename,
    #                       activators=activators)
    elif plot_type == 'release_endpoint_no_norm':
        nba.plot_release_endpoints(df, nbd_residues, normalized_to_wt=False,
                           last_n_pts=3, file_basename=file_basename,
                           activators=activators)
    elif plot_type == 'initial_rate_samples':
        nba.plot_initial_rate_samples(df, nbd_residues, timepoint_ix=4,
                                      file_basename=file_basename,
                                      activators=activators)
    elif plot_type == 'fret_endpoint':
        nba.plot_nbd_endpoints(df_pre, nbd_residues, datatype='FRET',
                               last_n_pts=10, file_basename=file_basename,
                               activators=activators)
    # Unknown plot type
    else:
        print "Unknown plot type: %s" % plot_type
        sys.exit(1)

