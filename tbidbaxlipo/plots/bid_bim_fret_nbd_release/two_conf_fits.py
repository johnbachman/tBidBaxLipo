from tbidbaxlipo.plots.nbd_bax_analysis import plot_2conf_fits, plot_3conf_fits
from tbidbaxlipo.plots.bid_bim_fret_nbd_release.preprocess_data \
        import df_pre, df, nbd_residues
from matplotlib import pyplot as plt

plt.ion()
nbd_residues=['184']
fret_fr = plot_2conf_fits(df_pre, nbd_residues, activator='Bim', dtype='FRET',
                replicates=(1, 2, 3))
#plt.show()
nbd_fr = plot_3conf_fits(df_pre, nbd_residues, activator='Bim', dtype='NBD',
                replicates=(1, 2, 3))
#plt.show()

fret_k1 = [fret_fr[i].param_dict['c0_to_c1_k'] for i in range(3)]
nbd_k1 = [nbd_fr[i].param_dict['c0_to_c1_k'] for i in range(3)]


print "FRET", fret_k1
print "NBD", nbd_k1
