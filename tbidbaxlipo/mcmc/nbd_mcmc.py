"""
Functions for fitting the mechanistic NBD insertion models to the NBD-Bax mutant
fluorescence data using MCMC.

.. todo:: Likelihood should be scaled by number of timepoints, otherwise it

overwhelms the prior for "no reason".

.. todo:: Ideally, would have a way of pre-equilibrating the system for the

just-Bax condition, and then perturb it with the addition of tBid.
"""

import bayessb
from pysb.integrate import odesolve
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import tbidbaxlipo.mcmc
from tbidbaxlipo.util.report import Report
from matplotlib.font_manager import FontProperties
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

model_names = ['ta', 'tai', 'taid', 'taidt', 'tair', 'taird', 'tairdt',
               'tad', 'tadt', 'tar', 'tard', 'tardt']

nbd_site_names = ['c3', 'c62']

class NBD_MCMC(tbidbaxlipo.mcmc.MCMC):
    """Fit mechanistic tBid/Bax models to NBD data.

    Initialize parent bayessb.MCMC and then set additional fields.

    Parameters
    ----------
    options : MCMCOpts
        Options for MCMC initialization.
    nbd_avgs : numpy.array
        data mean
    nbd_stds : numpy.array
        data SD
    nbd_sites : list of strings
        Sites from data to fit.
    nbd_observables : list of strings
        Observables from model to fit to the sites in nbd_sites.
    builder : tbidbaxlipo.models.core.Builder
    """

    def __init__(self, options, nbd_avgs, nbd_stds, nbd_sites, nbd_observables,
                 builder):
        # Call the superclass constructor
        tbidbaxlipo.mcmc.MCMC.__init__(self, options, builder)

        # Set the NBD-specific fields
        self.nbd_avgs = nbd_avgs
        self.nbd_stds = nbd_stds
        self.nbd_sites = nbd_sites
        self.nbd_observables = nbd_observables

        # Set the MCMC functions
        self.options.likelihood_fn = self.get_likelihood_function()
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    # NBD Functions
    # ==============
    def generate_figures(self, report_name='report', do_report=True,
                         mixed_start=None, num_samples=500):
        """Plots a series of useful visualizations of the walk, the
        quality of the fit, etc."""
        plt.ion()

        if do_report:
            rep = Report()

        # Plot "Before" curves -------
        before_fig = self.fit_plotting_function(self.initial_position)
        if do_report:
            rep.add_figure(before_fig)

        # Print initial parameter values
        init_vals = zip([p.name for p in self.options.model.parameters],
                          self.cur_params(position=self.initial_position))
        init_vals_str = 'Initial values:\n'
        init_vals_str += '\n'.join(['%s: %g' % (init_vals[i][0],
                                                init_vals[i][1])
                                     for i in range(0, len(init_vals))])
        print init_vals_str
        if do_report:
            rep.add_text(init_vals_str)

        # Plot "After" curves ------------
        # Set to last fit position
        if mixed_start is None:
            mixed_start = self.options.nsteps / 2
        mixed_positions = self.positions[mixed_start:,:]
        mixed_accepted_positions = mixed_positions[self.accepts[mixed_start:]]
        last_position = mixed_accepted_positions[-1,:]

        after_fig = self.fit_plotting_function(last_position)
        if do_report:
            rep.add_figure(after_fig)

        # Print final parameter values
        last_fit_params = self.cur_params(position=last_position)
        last_vals = zip([p.name for p in self.options.model.parameters],
                           last_fit_params)
        last_vals_str = 'Final values:\n'
        last_vals_str += '\n'.join(['%s: %g' % (last_vals[i][0],
                                                last_vals[i][1])
                                    for i in range(0, len(last_vals))])
        print last_vals_str
        if do_report:
            rep.add_text(last_vals_str)

        """
        # Plot sampling of fits # FIXME should be real sampling------
        plt.figure()
        plot_data(nbd_site)
        max_position_index = len(mixed_accepted_positions) - 1
        num_to_plot = min(num_samples, max_position_index)
        for i in range(num_to_plot):
            cur_position = mixed_accepted_positions[max_position_index - i,:]
            x = nbd_timecourse(mcmc, cur_position, nbd_observable)
            plt.plot(mcmc.options.tspan, x, 'b', label=nbd_site, alpha=0.02)
        plt.title("Last %d accepted positions" % num_to_plot)
        plt.show()
        if do_report:
            rep.add_current_figure()

        # Plot marginal distributions of parameters ---------------
        for i, cur_param in enumerate(mcmc.options.estimate_params):
            plt.figure()
            plt.hist(mixed_accepted_positions[:,i])
            plt.title("%s, last %d accepts: initval %f" %
                    (cur_param.name, len(mixed_accepted_positions[:,i]),
                     cur_param.value))
            plt.show()

            if do_report:
                rep.add_current_figure()

        # Plot convergence traces of all parameters
        plt.figure()
        plt.plot(mcmc.positions)
        plt.title("Parameter traces")
        plt.legend([p.name for p in mcmc.options.estimate_params],
                    loc='lower left', prop={'size':7})
        plt.show()
        if do_report:
            rep.add_current_figure()
        """

        # Add code to report
        if do_report:
            # Add the code for the fitting (this file)
            rep.add_python_code('nbd_mcmc_pysb.py')

            # Add the code for the model
            rep.add_python_code('models/core.py')
            rep.add_python_code('models/one_cpt.py')

            # Write the report 
            rep.write_report(report_name)

    def plot_data(self, axis):
        """Plots the current data into the given axis."""
        alpha = 0.5
        if 'c3' in self.nbd_sites:
            axis.plot(self.options.tspan, self.nbd_avgs[0], 'r.',
                     label='c3 data', alpha=alpha)
        if 'c62' in self.nbd_sites:
            axis.plot(self.options.tspan, self.nbd_avgs[1], 'g.',
                     label='c62 data', alpha=alpha)
        #plt.plot(nbd.time_other, nbd_avgs[2], 'b.', label='c120 data',
        #          alpha=alpha)
        #plt.plot(nbd.time_other, nbd_avgs[3], 'm.', label='c122 data',
        #         alpha=alpha)
        #plt.plot(nbd.time_other, nbd_avgs[4], 'k.', label='c126 data',
        #         alpha=alpha)

    def get_observable_timecourses(self, position):
        """Gets the timecourses associated with the experimental observables.

        Parameters
        ----------
        position : numpy.array
            Values of parameters (in log10) at desired position.

        Returns
        -------
        dict
            Dict containing the a timecourse for every observable. Keys are
            the observable names.
        """
        # TODO Need to be able to get the indices for scaling parameters
        # from the model so that they're not hardcoded

        # Run the simulation and scale the simulated timecourses
        yout = self.simulate(position, observables=True)
        params = self.cur_params(position)
        timecourses = {}
        for obs in self.nbd_observables:
            timecourses[obs] = ((yout[obs] /
                           self.options.model.parameters['Bax_0'].value)
                          * params[3])

        return timecourses

    # A function to generate the likelihood function
    def get_likelihood_function(self):
        """Returns a likelihood function for the specified NBD site."""
        data_indices = []
        if 'c3' in self.nbd_sites:
            data_indices.append(0)
        if 'c62' in self.nbd_sites:
            data_indices.append(1)
        # Make sure the list is not empty
        if not data_indices:
            raise Exception('Failed to initialize data_indices!')

        # Make sure that there is a corresponding data trajectory for each
        # observable
        if len(data_indices) != len(self.nbd_observables):
            raise Exception('Length of list of nbd_sites does not match the '
                            'list of nbd_observables!')

        def likelihood(mcmc, position):
            yout = mcmc.simulate(position, observables=True)
            # TODO Need to be able to get the indices from the model so that
            # they're not hardcoded
            params = mcmc.cur_params(position)
            err = 0
            for data_index, obs_name in zip(data_indices, self.nbd_observables):
                timecourse = ((yout[obs_name] /
                               mcmc.options.model.parameters['Bax_0'].value)
                              * params[3])
                err += np.sum((self.nbd_avgs[data_index] - timecourse)**2 /
                                 (2 * self.nbd_stds[data_index]**2))
            return err

        return likelihood

    def get_basename(self):
        """A function for standardizing, in one place, the format for pickled
        NBD_MCMC objects.
        """
        return '%s_%s_%s_%d_T%.2f_s%d' % (self.options.model.name,
                                    '-'.join(self.nbd_sites),
                                    '-'.join(self.nbd_observables),
                                    self.options.nsteps,
                                    self.options.T_init,
                                    #np.log10(self.options.thermo_temp),
                                    self.options.seed)


