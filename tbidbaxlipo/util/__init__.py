import matplotlib

colors = ['r', 'g', 'b', 'c', 'y', 'm', 'k']

def color_iter():
    return iter(colors)
# A colorblind safe, print friendly, and photocopy safe color scheme:
# http://colorbrewer2.org/?type=diverging&scheme=PuOr&n=3
fig_orange = '#f1a340'
fig_gray = '#f7f7f7'
fig_purple = '#998ec3'

def set_fig_params_for_publication():
    matplotlib.rcParams['font.size'] = 6
    matplotlib.rcParams['font.sans-serif'] = 'Arial'
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['text.latex.preamble'] = [
            r'\usepackage{helvet}',
            r'\usepackage{sansmath}',
            r'\sansmath']
    matplotlib.rcParams['xtick.major.size'] = 2
    matplotlib.rcParams['ytick.major.size'] = 2
    matplotlib.rcParams['xtick.major.pad'] = 2
    matplotlib.rcParams['ytick.major.pad'] = 2

