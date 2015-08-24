import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d

#draw a vector
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

from tbidbaxlipo.util import set_fig_params_for_publication

# Credit for the arrow code goes to
# http://stackoverflow.com/questions/11140163/
#           python-matplotlib-plotting-a-3d-cube-a-sphere-and-a-vector
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


def plot_slice_diagram(fixed_concs):
    """Plots a set of 3D axes with a super-imposed, to serve as a schematic
    for experimental designs where concentrations are varied across two axes
    while the third is held fixed.

    Parameters
    ----------
    fixed_concs : one- or two-tuple
        If there is one element in the tuple, then one axis is held fixed and a
        plane is plotted. If there are two elements, then two axes are held
        fixed and a line is plotted. The elements of the tuple should be either
        'x', 'y', or 'z', specifying which of the three axes is held fixed.
    """

    xx, yy = np.meshgrid(range(10), range(10))

    z = np.zeros((10, 10))
    z.fill(5)

    ax = plt.figure(figsize=(1.0, 1.0), dpi=300).gca(projection='3d')

    if len(fixed_concs) == 1:
        if fixed_concs[0] == 'x':
            ax.plot_surface(z, yy, xx, color='r', alpha=0.5)
        elif fixed_concs[0] == 'y':
            ax.plot_surface(xx, z, yy, color='r', alpha=0.5)
        elif fixed_concs[0] == 'z':
            ax.plot_surface(xx, yy, z, color='r', alpha=0.5)
        else:
            raise ValueError('The argument fixed_conc must be one '
                             'of ["x", "y", "z"]')
    elif len(fixed_concs) == 2:
        if 'x' in fixed_concs and 'y' in fixed_concs:
            a = Arrow3D([0.5, 0.5], [0.5, 0.5], [0, 1.], mutation_scale=10,
                        lw=1, arrowstyle="-|>", color='k')
        elif 'x' in fixed_concs and 'z' in fixed_concs:
            a = Arrow3D([0.5, 0.5], [0, 1.0], [0.5, 0.5], mutation_scale=10,
                        lw=1, arrowstyle="-|>", color='k')
        elif 'y' in fixed_concs and 'z' in fixed_concs:
            a = Arrow3D([0, 1.0], [0.5, 0.5], [0.5, 0.5], mutation_scale=10,
                        lw=1, arrowstyle="-|>", color='k')
        else:
            raise ValueError('The elements in fixed_concs must be one '
                             'of ["x", "y", "z"]')
        ax.add_artist(a)
    else:
        raise ValueError('fixed_concs must be of length 1 or 2.')

    #fig.plot_surface(xx, yy, zbottom, color='r')
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
    #plt.ion()

    plot_slice_diagram(fixed_concs=('x', 'y'))
    plt.savefig('slice_bid_bax_fixed.pdf')

    plot_slice_diagram(fixed_concs=('y', 'z'))
    plt.savefig('slice_bax_lipos_fixed.pdf')

    plot_slice_diagram(fixed_concs=('x', 'z'))
    plt.savefig('slice_bid_lipos_fixed.pdf')

    plot_slice_diagram(fixed_concs='x')
    plt.savefig('slice_bid_fixed.pdf')

    plot_slice_diagram(fixed_concs='y')
    plt.savefig('slice_bax_fixed.pdf')

    plot_slice_diagram(fixed_concs='z')
    plt.savefig('slice_lipos_fixed.pdf')

