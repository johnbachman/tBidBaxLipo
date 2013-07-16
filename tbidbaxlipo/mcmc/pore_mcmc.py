"""
A class for fitting the mechanistic pore formation models to dye release data.
"""

import tbidbaxlipo.mcmc
from tbidbaxlipo.models import lipo_sites, one_cpt
import bayessb
import pickle
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import collections
import sys

model_names = ['bax_heat',
               'bax_heat_reversible',
               'bax_heat_dimer',
               'bax_heat_dimer_reversible',
               'bax_heat_auto',
               'bax_heat_auto_reversible',
               'bax_heat_auto_dimer',
               'bax_heat_auto_dimer_reversible',
              ]

class PoreMCMC(tbidbaxlipo.mcmc.MCMC):
    """Fit mechanistic tBid/Bax models to dye release titration data.

    Initialize parent tbidbaxlipo.mcmc.MCMC and then set additional fields.

    Parameters
    ----------
    options : MCMCOpts
        Options for MCMC initialization.
    data :
        pandas data structure containing the timecourses for each Bax condition.
    builder : tbidbaxlipo.models.core.Builder
    """
    def __init__(self, options, data, dataset_name, builder):
        # Call the superclass constructor
        tbidbaxlipo.mcmc.MCMC.__init__(self, options, builder)

        # Store the data
        self.data = data
        self.dataset_name = dataset_name

        # Set the MCMC functions
        self.options.likelihood_fn = self.get_likelihood_function()
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    # Implementations of necessary functions
    def get_likelihood_function(self):

        def likelihood(mcmc, position):
            err = 0
            for bax_conc in mcmc.data.columns:
                # Get the data for this concentration
                tc = mcmc.data[bax_conc]
                y_data  = np.array(tc[:,'MEAN'])
                time = np.array(tc[:,'TIME'])
                mcmc.solver.tspan = time # set the time span

                # Get the simulated data for this concentration
                mcmc.options.model.parameters['Bax_0'].value = bax_conc
                x = mcmc.simulate(position=position, observables=True)
                avg_pores = x['pores']/ \
                            mcmc.options.model.parameters['Vesicles_0'].value
                y_mod = 1 - np.exp(-avg_pores)

                # Calculate the error, accounting for the SD at this
                # concentration.
                # Skip the first timepoint--the SD is 0 (due to normalization)
                # and hence gives nan when calculating the error.
                err += np.sum(((y_data[1:] - y_mod[1:])**2) / \
                        (2 * (np.array(tc[:,'SD'][1:]) ** 2)))
            return err
        return likelihood

    def plot_data(self, axis):
        # Plot the titration of Bax timecourses
        for bax_conc in self.data.columns:
            tc = self.data[bax_conc]
            axis.plot(tc[:,'TIME'], tc[:,'MEAN'], # error=tc[:,'SD'],
                       color='gray')

    def get_observable_timecourses(self, position):
        """Return the timecourses for all concentrations."""
        timecourses = collections.OrderedDict()

        for bax_conc in self.data.columns:
            # Get the timepoints for this concentration
            tc = self.data[bax_conc]
            time = np.array(tc[:,'TIME'])
            self.solver.tspan = time # set the time span
            self.options.model.parameters['Bax_0'].value = bax_conc
            x = self.simulate(position=position, observables=True)
            avg_pores = x['pores'] / \
                        self.options.model.parameters['Vesicles_0'].value
            y_mod = 1 - np.exp(-avg_pores)
            timecourses['Bax %d nM' % bax_conc] = [time, y_mod]
        return timecourses

    def get_basename(self):
        return '%s_%s_%s_%d_s%d' % (self.dataset_name,
                                 self.builder.get_module(),
                                 self.options.model.name,
                                 self.options.nsteps,
                                 self.options.seed)

def get_PoreMCMC_instance():
    b = one_cpt.Builder()
    b.build_model_bax_heat()

    data = pickle.load(open('../data/130614_norm_timecourses_df.pck'))
    dataset_name = '130614'

    opts = bayessb.MCMCOpts()
    opts.model = b.model
    opts.tspan = np.linspace(0, 6000, 196)
    opts.estimate_params = [p for p in b.model.parameters
                            if not p.name.endswith('_0')]
    opts.initial_values = [p.value for p in opts.estimate_params]
    opts.nsteps = 10
    opts.T_init = 1
    opts.anneal_length = 0
    opts.use_hessian = False
    opts.sigma_step = 0
    opts.norm_step_size = 0.1
    opts.seed = 1

    pm = PoreMCMC(opts, data, dataset_name, b)
    return pm

def test_PoreMCMC_init():
    pm = get_PoreMCMC_instance()
    assert True

def test_plot_data():
    pm = get_PoreMCMC_instance()
    fig = plt.figure()
    axis = fig.gca()
    pm.plot_data(axis)
    #plt.show()
    assert True

def test_plot_fit():
    pm = get_PoreMCMC_instance()
    pm.initialize()
    fig = pm.fit_plotting_function(position=pm.initial_position)
    fig.savefig('test_plot_fit.png')

if __name__ == '__main__':
    # Allowable compartmentalization type names:
    cpt_types = ['one_cpt', 'lipo_sites']

    # Prepare the data
    # ================
    data = pickle.load(open('../data/130614_norm_timecourses_df.pck'))
    dataset_name = '130614'

    # Parse the args
    # ==============
    # Keyword args are set at the command line as e.g., key=val
    # and subsequently split at the equals sign
    kwargs = dict([arg.split('=') for arg in sys.argv[1:]])

    # We set these all to None so later on we can make sure they were
    # properly initialized.
    random_seed = None
    model = None

    print "Keyword arguments: "
    print kwargs

    # Before we begin, we make sure we have all the keyword arguments that
    # we are going to need.
    if 'random_seed' not in kwargs or \
       'cpt_type' not in kwargs or \
       'nsteps' not in kwargs or \
       'model' not in kwargs:
        raise Exception('One or more needed arguments was not specified! ' \
                'Arguments must include random_seed, model, cpt_type, ' \
                'and nsteps.')

    # First we get the compartmentalization type:
    cpt_type = kwargs['cpt_type']
    if cpt_type not in cpt_types:
        raise Exception('Allowable values for cpt_type are: one_cpt, '
                        'lipo_sites')

    # Get the appropriate builder
    if cpt_type == 'lipo_sites':
        builder = lipo_sites.Builder()
    elif cpt_type == 'one_cpt':
        builder = one_cpt.Builder()

    # Now we get the model, which is specified as a string from the
    # set seen below.
    if kwargs['model'] not in model_names:
        raise Exception('%s is not an allowed model!' %
                        kwargs['model'])
    else:
        # Here we use a bit of Python trickery to avoid a long list of
        # if/elif statements: we append the model abbreviation to the
        # build_model function name and then eval it:
        build_model_cmd = 'builder.build_model_%s()' % kwargs['model']
        eval(build_model_cmd)
        model = builder.model

    # Set the random seed:
    random_seed = int(kwargs['random_seed'])

    # Set the number of steps:
    nsteps = int(kwargs['nsteps'])

    # Set the initial temperature
    if 'T_init' in kwargs:
        T_init = float(kwargs['T_init'])
    else:
        T_init = 1

    # A sanity check to make sure everything worked:
    if None in [model, random_seed, nsteps]:
        raise Exception('Something went wrong! One of the arguments to ' \
                        'do_fit was not initialized properly.')

    # We set the random_seed here because it affects our choice of initial
    # values
    np.random.seed(random_seed)

    # Initialize the MCMC arguments
    opts = bayessb.MCMCOpts()
    opts.model = model
    opts.tspan = np.linspace(0, 6000, 196) # This should not be used
    opts.estimate_params = builder.estimate_params
    opts.initial_values = builder.random_initial_values()
    opts.nsteps = nsteps
    opts.norm_step_size = 0.01
    opts.sigma_step = 0
    #opts.sigma_max = 50
    #opts.sigma_min = 0.01
    #opts.accept_rate_target = 0.23
    #opts.accept_window = 100
    #opts.sigma_adj_interval = 200
    opts.anneal_length = 0
    opts.use_hessian = False #True # TRUE
    opts.hessian_scale = 1
    opts.hessian_period = opts.nsteps / 20 #10
    opts.seed = random_seed
    opts.T_init = 1 # T_init
    mcmc = PoreMCMC(opts, data, dataset_name, builder)

    mcmc.do_fit()

    # Pickle it
    output_file = open('%s.mcmc' % mcmc.get_basename(), 'w')
    pickle.dump(mcmc, output_file)
    output_file.close()

    print "Done."

