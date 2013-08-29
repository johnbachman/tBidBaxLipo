from bayessb.report import reporter, Result, ThumbnailResult
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from tbidbaxlipo.util import fitting
from tbidbaxlipo.plots.bax_heat import fit_from_dataframe, two_exp_func

reporter_group_name = 'Bax activated by heat'
num_samples = 10

two_exp_fits_dict = {}

@reporter('Two-exponential fits')
def two_exponential_fits(mcmc_set):
    # Use the first MCMC in the chain to get the data, etc.
    mcmc = mcmc_set.chains[0]

    # Run a simulation
    (max_lkl, max_lkl_position) = mcmc_set.maximum_likelihood()
    data = mcmc.get_observables_as_dataframe(max_lkl_position)
    bax_concs = np.array(data.columns, dtype='float')

    plot_filename = '%s_two_exp_fits.png' % (mcmc_set.name)
    thumbnail_filename = '%s_two_exp_fits_th.png' % (mcmc_set.name)

    (fmax_arr, k1_arr, k2_arr) = fit_from_dataframe(data)

    fig = Figure()
    for i, bax_conc in enumerate(bax_concs):
        ax = fig.gca()
        conc_data = data[bax_conc]
        time = conc_data[:, 'TIME']
        y = conc_data[:, 'MEAN']

        two_exp_fit = two_exp_func(time, fmax_arr[i], k1_arr[i], k2_arr[i])
        # Plot original sim
        ax.plot(time, y, 'b')
        # Plot fitted func
        ax.plot(time, two_exp_fit, 'r')

    ax.set_title('Two-exponential fits for %s' % mcmc_set.name)
    ax.set_xlabel('Time')
    ax.set_ylabel('Dye release')
    #ax.legend(loc='lower right')
    canvas = FigureCanvasAgg(fig)
    fig.set_canvas(canvas)
    fig.savefig(plot_filename)
    fig.savefig(thumbnail_filename, dpi=10)

    two_exp_fits_dict[mcmc_set.name] = (fmax_arr, k1_arr, k2_arr)
    return ThumbnailResult(thumbnail_filename, plot_filename)

@reporter('Fmax curve')
def fmax_curve(mcmc_set):
    return plot_parameter_curve(mcmc_set, 0, 'Fmax')

@reporter('k1 curve')
def k1_curve(mcmc_set):
    return plot_parameter_curve(mcmc_set, 1, 'k1')

@reporter('k2 curve')
def k2_curve(mcmc_set):
    return plot_parameter_curve(mcmc_set, 2, 'k2')

def plot_parameter_curve(mcmc_set, p_index, p_name):
    # Make sure we've already run the fits for this mcmc set!
    if mcmc_set.name not in two_exp_fits_dict.keys():
        raise Exception('%s not found in two_exp_fits_dict!' % mcmc_set.name)
    # Make sure we've already run the fits for the data!
    if 'data' not in two_exp_fits_dict.keys():
        fit_data(mcmc_set)

    # Get the parameter array
    p_arr = two_exp_fits_dict[mcmc_set.name][p_index]
    p_arr_data = two_exp_fits_dict['data'][p_index]
    data = mcmc_set.chains[0].data
    plot_filename = '%s_%s_curve.png' % (mcmc_set.name, p_name)
    thumbnail_filename = '%s_%s_curve_th.png' % (mcmc_set.name, p_name)

    # Plot of parameter
    fig = Figure()
    ax = fig.gca()
    ax.plot(data.columns, p_arr, 'b')
    ax.plot(data.columns, p_arr_data, marker='o', linestyle='', color='r')
    ax.set_ylabel('%s value' % p_name)
    ax.set_xlabel('[Bax] (nM)')
    ax.set_title('%s for %s' % (mcmc_set.name, p_name))
    canvas = FigureCanvasAgg(fig)
    fig.set_canvas(canvas)
    fig.savefig(plot_filename)
    fig.savefig(thumbnail_filename, dpi=10)
    return ThumbnailResult(thumbnail_filename, plot_filename)

def fit_data(mcmc_set):
    data = mcmc_set.chains[0].data
    fmax_arr = np.zeros(len(data.columns))
    k1_arr = np.zeros(len(data.columns))
    k2_arr = np.zeros(len(data.columns))
    bax_concs = np.array(data.columns, dtype='float')
    for i, bax_conc in enumerate(bax_concs):
        conc_data = data[bax_conc]
        time = conc_data[:, 'TIME']
        y = conc_data[:, 'MEAN']
        fmax = fitting.Parameter(0.9)
        k1 = fitting.Parameter(0.01)
        k2 = fitting.Parameter(0.001)
        def two_part_exp(t):
            return (fmax() * (1 - np.exp(-k1() * (1 - np.exp(-k2()*t)) * t)))
        fitting.fit(two_part_exp, [fmax, k1, k2], y, time)
        fmax_arr[i] = fmax()
        k1_arr[i] = k1()
        k2_arr[i] = k2()
    two_exp_fits_dict['data'] = (fmax_arr, k1_arr, k2_arr)


