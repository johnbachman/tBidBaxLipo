from tbidbaxlipo.plots.nbd_bax_analysis import plot_2conf_fits, plot_3conf_fits
from tbidbaxlipo.plots.x141203_Bax_Bax_FRET.preprocess_data \
        import df_pre, df, nbd_residues
from matplotlib import pyplot as plt
plt.ion()

plot_2conf_fits(df_pre, nbd_residues, activator='Bid', dtype='FRET',
                replicates=(1,))

#plot_3conf_fits(df_pre, nbd_residues, activator='Bid', dtype='FRET',
#                replicates=(1,))

