from matplotlib import pyplot as plt

line_colors = {'Bid': 'r', 'Bim': 'b'}
line_styles = {1:':', 2:'-', 3:'--'}
rep_colors = {1:'r', 2:'g', 3:'b'}

def plot_all(df, nbd_residues, file_basename=None):
    for nbd_index, nbd_site in enumerate(nbd_residues):
        plt.figure(figsize=(11, 5))
        # Make the release plot
        plt.subplot(1, 2, 1)
        activators = ['Bid', 'Bim']
        for activator in activators:
            for i in range(1, 4):
                t = df[(activator, 'Release', nbd_site, i, 'TIME')]
                v = df[(activator, 'Release', nbd_site, i, 'VALUE')]
                plt.plot(t, v, label='%s Rep %d' % (activator, i),
                        color=line_colors[activator],
                        linestyle=line_styles[i])

                plt.xlabel('Time (sec)')
                plt.ylabel('Pct. Release')
                plt.title('Release for NBD-%s-Bax' % nbd_site)
                plt.legend(loc='lower right')
        # There is no NBD curve for WT Bax, so skip the NBD
        # plot
        if nbd_site == 'WT':
            continue
        # Make the NBD plot
        plt.subplot(1, 2, 2)
        for activator in ['Bid', 'Bim']:
            for i in range(1, 4):
                t = df[(activator, 'NBD', nbd_site, i, 'TIME')]
                v = df[(activator, 'NBD', nbd_site, i, 'VALUE')]
                plt.plot(t, v, label='%s Rep %d' % (activator, i),
                        color=line_colors[activator],
                        linestyle=line_styles[i])
                plt.xlabel('Time (sec)')
                plt.ylabel('F/F0')
                plt.title('F/F0 for NBD-%s-Bax' % nbd_site)
                plt.legend(loc='lower right')
        plt.tight_layout()
        if file_basename:
            plt.savefig('%s_%s.pdf' % (file_basename, nbd_index))
            plt.savefig('%s_%s.png' % (file_basename, nbd_index))

if __name__ == '__main__':
    from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
    plot_all(df, nbd_residues)
