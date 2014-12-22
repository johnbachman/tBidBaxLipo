from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
from tbidbaxlipo.util import set_fig_params_for_publication

def plot_slice_diagram(fixed_conc):
    """Plots a set of 3D axes with a super-imposed, to serve as a schematic
    for experimental designs where concentrations are varied across two axes
    while the third is held fixed.

    Parameters
    ----------
    fixed_conc : string
        One of 'x', 'y', or 'z', specifying which of the three axes is held
        fixed.
    """

    xx, yy = np.meshgrid(range(10), range(10))

    z = np.zeros((10, 10))
    z.fill(5)

    plt.ion()
    ax = plt.figure(figsize=(1.0, 1.0), dpi=300).gca(projection='3d')
    if fixed_conc == 'x':
        ax.plot_surface(z, yy, xx, color='r', alpha=0.5)
    elif fixed_conc == 'y':
        ax.plot_surface(xx, z, yy, color='r', alpha=0.5)
    elif fixed_conc == 'z':
        ax.plot_surface(xx, yy, z, color='r', alpha=0.5)
    else:
        raise ValueError('The argument fixed_conc must be one '
                         'of ["x", "y", "z"]')

    #fig.plot_surface(xx, yy, zbottom, color='r')
    ax = plt.gca()
    ax.set_xlabel('[cBid]')
    ax.set_ylabel('[Bax]')
    ax.set_zlabel('[Liposomes]')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_frame_on(False)
    ax.tick_params(which='both', pad=0, length=0)
    ax.xaxis.labelpad = 0
    ax.xaxis.label.set_size(6)
    ax.yaxis.label.set_size(6)
    ax.zaxis.label.set_size(6)
    ax.view_init(elev=25., azim=-130.)
    ax.xaxis._axinfo['label']['space_factor'] = 0.9
    ax.yaxis._axinfo['label']['space_factor'] = 0.9
    ax.zaxis._axinfo['label']['space_factor'] = 0.9
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1)

if __name__ == '__main__':

    plot_slice_diagram(fixed_conc='x')
    plt.savefig('slice_bid_fixed.pdf')
    plot_slice_diagram(fixed_conc='y')
    plt.savefig('slice_bax_fixed.pdf')
    plot_slice_diagram(fixed_conc='z')
    plt.savefig('slice_lipos_fixed.pdf')

