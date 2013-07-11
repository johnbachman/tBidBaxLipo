from matplotlib import pyplot as plt

def plot_raw(data, display=False):
    # For each mutant, make a plot of all replicates
    nbd_names = data.columns.levels[0]
    replicates = data.columns.levels[1]

    if display:
        plt.ion()

    for nbd_name in nbd_names:
        plt.figure()
        for replicate in replicates:
            tc = data[(nbd_name, replicate)]
            plt.plot(tc[:, 'TIME'], tc[:, 'VALUE'],
                     label='%s rep. %d' % (nbd_name, replicate))
            plt.legend(loc='lower right')
            plt.xlabel('Time (sec)')
            plt.ylabel('Fluorescence (AU)')
            plt.title(nbd_name)
        if display:
            plt.show()

    if display:
        plt.ioff()

