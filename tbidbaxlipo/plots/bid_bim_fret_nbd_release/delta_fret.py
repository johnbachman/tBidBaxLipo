from tbidbaxlipo.plots import nbd_bax_analysis as nba

from tbidbaxlipo.plots.bid_bim_fret_nbd_release.preprocess_data \
            import df_pre, nbd_residues
from matplotlib import pyplot as plt
from tbidbaxlipo.models.nbd.multiconf import Builder
import numpy as np
import cPickle

plt.ion()

# Dict to store the results
delta_dict = {}

residues = [res for res in nbd_residues if res != '62']

# --- First, fit Bid data ---
"""
bid_fit_results = nba.plot_3conf_fits(df_pre, residues, 'Bid', dtype='FRET',
                                      do_plot=False)

for fr in bid_fit_results:
    # Calculate the difference
    key = (fr.activator, fr.measurement, fr.nbd_site, fr.rep_index)
    delta_fret = np.max(fr.y) - fr.y[-1]
    delta_dict[key] = delta_fret

with open('fret_deltas_bid.pck', 'w') as f:
    cPickle.dump(delta_dict, f)
"""

# --- Now do Bim ---
bim_fit_results = nba.plot_3conf_fits(df_pre, residues, 'Bim', dtype='FRET',
                                      do_plot=False)

for fr in bim_fit_results:
    # Calculate the difference
    key = (fr.activator, fr.measurement, fr.nbd_site, fr.rep_index)
    delta_fret = np.max(fr.y) - fr.y[-1]
    delta_dict[key] = delta_fret

with open('fret_deltas_bim.pck', 'w') as f:
    cPickle.dump(delta_dict, f)


