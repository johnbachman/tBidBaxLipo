import numpy as np
from scipy import optimize
from pysb.integrate import Solver
from matplotlib import pyplot as plt
import emcee
from emcee.utils import MPIPool
import mpi4py
import sys
from scipy.stats import pearsonr
import cPickle
import warnings
from tbidbaxlipo.models import one_cpt, lipo_sites
from tbidbaxlipo.models.nbd import multiconf

def posterior(position, gf):
    """A generic log posterior function."""
    post_val = prior(position, gf) + likelihood(position, gf)
    return post_val

def negative_posterior(position, gf):
    """A generic negative log posterior function.

    Use the negative log posterior when using a minimization (rather than
    probability maximization) algorithm.
    """
    post_val = prior(position, gf) + likelihood(position, gf)
    print("posterior: %s" % post_val)
    return -post_val

def prior(position, gf):
    """A generic prior function."""
    prior_prob = 0
    for i, prior in enumerate(gf.priors):
        prior_prob += prior.pdf(position[i])
    return -prior_prob

def likelihood(position, gf):
    # A generic objective function
    # The cumulative error over all timecourses
    err = 0

    # Iterate over each condition (1st dimension of the data matrix)
    for cond_ix in range(gf.data.shape[0]):
        # Create a placeholder for a time offset, if there is one
        timeoffset = None
        # Set the parameters appropriately for the simulation:
        # Iterate over the globally fit parameters
        for g_ix, p in enumerate(gf.builder.global_params):
            p.value = 10 ** position[g_ix]
            if p.name == 'timeoffset':
                timeoffset = 10 ** position[g_ix]
        # Iterate over the locally fit parameter_s
        for l_ix, p in enumerate(gf.builder.local_params):
            ix_offset = len(gf.builder.global_params) + \
                        cond_ix * len(gf.builder.local_params)
            p.value = 10 ** position[l_ix + ix_offset]
        # Now fill in the initial condition parameters
        if gf.params is not None:
            for p_name, values in gf.params.iteritems():
                p = gf.builder.model.parameters[p_name]
                p.value = values[cond_ix]
        # Reset the timespan by adding one additional pt at the beginning
        if timeoffset:
            gf.solver.tspan = np.insert(gf.time, 0, -timeoffset)
        # Now run the simulation
        gf.solver.run()
        # Calculate the squared error over all the observables
        for obs_ix, obs_name in enumerate(gf.obs_name):
            if gf.use_expr:
                ysim = gf.solver.yexpr[obs_name]
            else:
                ysim = gf.solver.yobs[obs_name]
            # If we're using a time offset, skip the first point (the offset)
            # for the purposes of comparing to data
            if timeoffset:
                ysim = ysim[1:]
            # If integrator fails to converge, the results will contain NaN
            if np.any(np.isnan(ysim)):
                err = -np.inf
                continue
            # Get the data slice we want
            data = gf.data[cond_ix, obs_ix, :]
            # Get the appropriate SD for this data slice
            sigma = gf.data_sigma[cond_ix, obs_ix]
            # Calculate the log-likelihood
            loglkl = ((data - ysim) ** 2) / (2. * sigma ** 2)
            # Filter out the NaNs...
            filt_loglkl = loglkl[~np.isnan(loglkl)]
            # Take the sum
            err += -np.sum(filt_loglkl)

    return err

class GlobalFit(object):
    """Fit of PySB model to a set of multiple timecourses, with a
    mix of globally and locally fit parameters.

    Parameters
    ----------
    builder : pysb.builder.Builder
        Builder containing the model to fit. Should contain an attribute
        builder.global_params for the parameters that are to be fit globally.
    time : np.array
        The time vector.
    data : Three-dimensional np.array
        The experimental timecourses to fit.  The first dimension corresponds
        to the number of experimental conditions; the second dimension
        corresponds to the number of observables (in a given timecourse set);
        the third dimension corresponds to the timepoints.
    data_sigma : np.array
        Array of values with dimension corresponding to data indicating the
        standard deviation of the data.
    params : dict of lists, or None
        The keys to the dict should be names of parameters in the PySB model
        (e.g., initial conditions); each value should be a list containing
        values for the parameter for each of the entries in the data list. The
        length of each value in the dict should match the length of the data
        list. If None, indicates that there are no local initial conditions.
    obs_name : string
        The name of the model observable to compare against the data.
    obs_type : string, "Expression" or "Observable"
        Indicates whether the named expression/observable specified by
        obs_name is to be found in the model's set of Expression objects
        or Observable objects.

    Attributes
    ----------
    result : None or scipy.optimize.minimize fit result object
        The result field is initialized to None and is assigned the results
        of fitting after the :py:meth:`fit` method completes successfully.
    use_expr : boolean
        Based on the obs_type argument. True if the named observable is an
        Expression, False if it is an Observable.
    priors : list of priors
    solver : pysb.integrate.Solver
        A solver object used to run the model.
    """

    def __init__(self, builder, time, data, data_sigma, params, obs_name,
                 obs_type='Expression'):
        # Check that the dimensions of everything that has been provided matches
        # Check that the time vector matches the 3rd dimension of the data
        # vector
        if len(time) != data.shape[2]:
            raise ValueError("Length of time vector must match the length "
                             "of each data vector.")
        # Check that we don't have more than one condition but only set of
        # initial conditions
        if params is None and data.shape[0] != 1:
            raise ValueError("There are no initial condition parameters but "
                             "there is more than one condition in the data "
                             "matrix.")
        # Check that the number of initial conditions specified in the params
        # dict matches the first dimension of the data matrix
        if params is not None:
            for p, vals in params.iteritems():
                if not len(vals) == data.shape[0]:
                    raise ValueError("Each parameter in the params dict must "
                                     "have an entry for each entry in the "
                                     "data list.")
        # Check that the number of observables matches the 2nd dimension of
        # the data matrix
        if len(obs_name) != data.shape[1]:
            raise ValueError("The number of observables (%s) must match the "
                             "second dimension of the data matrix (%s)" %
                             (len(obs_name), data.shape[1]))
        # Check that there is a sigma in the data_sigma matrix for every
        # timecourse in data
        if data_sigma.shape[0] != data.shape[0] and \
           data_sigma.shape[1] != data.shape[1]:
            raise ValueError("data_sigma must specify an error SD for every "
                             "timecourse in the data matrix.")

        self.builder = builder
        self.time = time
        self.data = data
        self.data_sigma = data_sigma
        self.params = params
        self.obs_name = obs_name
        self.result = None
        if obs_type == 'Expression':
            self.use_expr = True
        elif obs_type == 'Observable':
            self.use_expr = False
        else:
            raise ValueError('obs_type must be Expression or Observable.')

        if self.builder.model.parameters.get('timeoffset'):
            use_time_offset = True
        else:
            use_time_offset = False

        self.init_solver(use_time_offset=use_time_offset)

        # Used to keep track of the number of steps run
        self.nstep = 0

        # Build up a list of priors corresponding to the global and local
        # parameters
        self.priors = []
        # Iterate over the globally fit parameters
        for g_ix, p in enumerate(self.builder.global_params):
            try:
                prior_index = self.builder.estimate_params.index(p)
                self.priors.append(self.builder.priors[prior_index])
            except ValueError:
                raise ValueError(
                        'The parameter %s, in global_params, must also be '
                        'present in estimate_params.' % p.name)
        # Iterate over the locally fit parameters
        for data_ix, data in enumerate(self.data):
            for l_ix, p in enumerate(self.builder.local_params):
                try:
                    prior_index = self.builder.estimate_params.index(p)
                    self.priors.append(self.builder.priors[prior_index])
                except ValueError:
                    raise ValueError(
                            'The parameter %s, in local_params, must also be '
                            'present in estimate_params.')

    def __getstate__(self):
        # Clear solver since it causes problems with pickling
        state = self.__dict__.copy()
        if 'solver' in state:
            del state['solver']
        return state

    def __setstate__(self, state):
        # Re-init the solver which we didn't pickle
        self.__dict__.update(state)
        if self.builder.model.parameters.get('timeoffset'):
            use_time_offset = True
        else:
            use_time_offset = False
        self.init_solver(use_time_offset=use_time_offset)

    def init_solver(self, use_time_offset=False):
        """Initialize solver from model and tspan."""
        Solver._use_inline = True
        # If we're using a time offset, note that it doesn't matter what value
        # goes in here, since it will be filled in by the fitting.
        if use_time_offset:
            tspan = np.insert(self.time, 0, 0)
        else:
            tspan = self.time
        self.solver = Solver(self.builder.model, tspan)

    def plot_func_single(self, x, data_ix, ax=None, alpha=1.0):
        x = 10 ** x

        s = Solver(self.builder.model, self.time)
        # Set the parameters appropriately for the simulation:
        # Iterate over the globally fit parameters
        for g_ix, p in enumerate(self.builder.global_params):
            p.value = x[g_ix]
        # Iterate over the locally fit parameters
        for l_ix, p in enumerate(self.builder.local_params):
            ix_offset = len(self.builder.global_params) + \
                        data_ix * len(self.builder.local_params)
            p.value = x[l_ix + ix_offset]
        # Now fill in the initial condition parameters
        for p_name, values in self.params.iteritems():
            p = self.builder.model.parameters[p_name]
            p.value = values[data_ix]
        # Now run the simulation
        s.run()
        # Plot the observable
        if ax is None:
            ax = plt.gca()
        if self.use_expr:
            ax.plot(self.time, s.yexpr[self.obs_name], color='r',
                     alpha=alpha)
        else:
            ax.plot(self.time, s.yobs[self.obs_name], color='r',
                     alpha=alpha)

    def plot_func(self, x, ax=None, obs_ix=0, plot_args=None,
                  normalize_to_f0=False):
        """Plots the timecourses with the parameter values given by x.

        Parameters
        ----------
        x : np.array or list
            The parameters to use for the plot, in log10 space.These should be
            in the same order used by the objective function: globally fit
            parameters first, then a set of local parameters for each of the
            timecourses being fit.
        normalize_to_f0: boolean
            If True, divides the observable trajectory by the initial value
            before plotting. Defaults to False.
        """
        if plot_args is None:
            plot_args = {}
        if ax is None:
            ax = plt.gca()
        # Iterate over each entry in the data array
        for cond_ix in range(self.data.shape[0]):
            self.set_parameters(x, obs_ix=obs_ix, cond_ix=cond_ix)
            # Now run the simulation
            self.solver.run()
            # Plot the observable
            obs_colors = ['r', 'g', 'b', 'k']
            obs_name = self.obs_name[obs_ix]
            # Get from either yexpr or yobs
            if self.use_expr:
                yplot = self.solver.yexpr[obs_name]
            else:
                yplot = self.solver.yobs[obs_name]
            # Normalize the curve if necessary
            if normalize_to_f0:
                yplot = yplot / float(yplot[0])
            # Plot
            ax.plot(self.solver.tspan, yplot, **plot_args)

    def set_parameters(self, x, obs_ix=0, cond_ix=0, plot_args=None):
        """Sets the parameter values in the model for simulation.

        Parameters
        ----------
        x : np.array or list
            The parameters to use, in log10 space.These should be
            in the same order used by the objective function: globally fit
            parameters first, then a set of local parameters for each of the
            timecourses being fit.
        """
        x = 10 ** x
        timeoffset = None
        # Set the parameters appropriately for the simulation:
        # Iterate over the globally fit parameters
        for g_ix, p in enumerate(self.builder.global_params):
            p.value = x[g_ix]
            if p.name == 'timeoffset':
                timeoffset = x[g_ix]
        # Iterate over the locally fit parameters
        for l_ix, p in enumerate(self.builder.local_params):
            ix_offset = len(self.builder.global_params) + \
                        cond_ix * len(self.builder.local_params)
            p.value = x[l_ix + ix_offset]
        # Now fill in the initial condition parameters
        if self.params is not None:
            for p_name, values in self.params.iteritems():
                p = self.builder.model.parameters[p_name]
                p.value = values[cond_ix]
        # Fill in the time offset, if there is one
        if timeoffset:
            self.solver.tspan = np.insert(self.time, 0, -timeoffset)

    def get_residuals(self, x, obs_ix=0, plot_args=None):
        """Gets the residuals with the parameter values given by x.

        Parameters
        ----------
        x : np.array or list
            The parameters to use for the simulation, in log10 space.
            These should be in the same order used by the objective function:
            globally fit parameters first, then a set of local parameters for
            each of the timecourses being fit.
        """
        if plot_args is None:
            plot_args = {}
        num_conditions = self.data.shape[0]
        num_timepoints = self.data.shape[2]
        residuals = np.zeros((num_conditions, num_timepoints))
        for cond_ix in range(num_conditions):
            self.set_parameters(x, obs_ix=obs_ix, cond_ix=cond_ix)
            # Run the simulation
            self.solver.run()
            # Get the simulated timecourse
            obs_name = self.obs_name[obs_ix]
            if self.use_expr:
                y = self.solver.yexpr[obs_name]
            else:
                y = self.solver.yobs[obs_name]
            # Calculate the residuals
            data = self.data[cond_ix, obs_ix, :]
            residuals[cond_ix, :] = data - y
        return residuals

def ens_sample(gf, nwalkers, burn_steps, sample_steps, threads=1,
               pos=None, random_state=None):
    """Samples from the posterior function using emcee.EnsembleSampler.

    The EnsembleSampler containing the chain is stored in gf.sampler.

    Note that parameters are log10-transformed during fitting, so the
    parameter values returned from the walk must be exponentiated to get
    them back to untransformed values (e.g., 10 ** gf.sampler.flatchain)

    Parameters
    ----------
    gf : emcee_fit.GlobalFit
        GlobalFit object containing the timepoints, data, builder object,
        Solver, etc.
    nwalkers : int
        Number of walkers to use in the emcee sampler.
    burn_steps : int
        Number of burn-in steps.
    sample_steps : int
        Number of sampling steps.
    threads : int
        Number of threads to use for parallelization. Default is 1.
    pos : numpy.array
        Matrix of initial positions for the chain. If None (default) random
        positions are chosen from the prior. Assigning a position allows
        previously run chains to be extended.
    random_state : random state for Mersenne Twister PRNG
        The random state to use to initialize the sampler's pseudo-random
        number generator. Can be used to continue runs from previous ones.
    """

    # Initialize the parameter array with initial values (in log10 units)
    # Number of parameters to estimate
    ndim = (len(gf.builder.global_params) +
            (len(gf.data) * len(gf.builder.local_params)))
    # Initialize the walkers with starting positions drawn from the priors
    # Note that the priors are in log10 scale already, so they don't
    # need to be transformed here
    if pos is None:
        p0 = np.zeros((nwalkers, ndim))
        for walk_ix in range(nwalkers):
            for p_ix in range(ndim):
                p0[walk_ix, p_ix] = gf.priors[p_ix].random()
    else:
        p0 = pos

    # Create the sampler object
    sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior,
                                         args=[gf],
                                         threads=threads)
    if random_state is not None:
        sampler.random_state = random_state

    print "Burn in sampling..."
    pos, prob, state = sampler.run_mcmc(p0, burn_steps, storechain=False)
    sampler.reset()

    print "Main sampling..."
    sampler.run_mcmc(pos, sample_steps)

    print "Done sampling."
    return sampler

def ens_mpi_sample(gf, nwalkers, burn_steps, sample_steps, pos=None,
                   random_state=None):
    pool = MPIPool(loadbalance=True)
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    # Initialize the parameter array with initial values (in log10 units)
    # Number of parameters to estimate
    ndim = (len(gf.builder.global_params) +
            (len(gf.data) * len(gf.builder.local_params)))
    # Initialize the walkers with starting positions drawn from the priors
    # Note that the priors are in log10 scale already, so they don't
    # need to be transformed here
    if pos is None:
        p0 = np.zeros((nwalkers, ndim))
        for walk_ix in range(nwalkers):
            for p_ix in range(ndim):
                p0[walk_ix, p_ix] = gf.priors[p_ix].random()
    else:
        p0 = pos

    # Create the sampler object
    sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior,
                                         args=[gf], pool=pool)
    if random_state is not None:
        sampler.random_state = random_state

    print "Burn in sampling..."
    pos, prob, state = sampler.run_mcmc(p0, burn_steps, storechain=False)
    sampler.reset()

    print "Main sampling..."
    sampler.run_mcmc(pos, sample_steps)

    # Close the pool!
    pool.close()

    print "Done sampling."
    return sampler

def pt_mpi_sample(gf, ntemps, nwalkers, burn_steps, sample_steps, thin=1,
                  pool=None, betas=None, pos=None, random_state=None,
                  pos_filename=None, convergence_interval=50):
    pool = MPIPool(loadbalance=True)
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    return pt_sample(gf, ntemps, nwalkers, burn_steps, sample_steps,
                     thin=thin, pool=pool, betas=betas, pos=pos,
                     random_state=random_state, pos_filename=pos_filename,
                     convergence_interval=convergence_interval)

def pt_sample(gf, ntemps, nwalkers, burn_steps, sample_steps, thin=1,
              pool=None, betas=None, pos=None, random_state=None,
              pos_filename=None, convergence_interval=50):
    """Samples from the posterior function.

    The emcee sampler containing the chain is stored in gf.sampler.

    Note that parameters are log10-transformed during fitting, so the
    parameter values returned from the walk must be exponentiated to get
    them back to untransformed values (e.g., 10 ** gf.sampler.flatchain)

    Parameters
    ----------
    ntemps : int
        The number of temperature to use in the temperature ladder.
    nwalkers : int
        Number of walkers to use in the emcee sampler.
    burn_steps : int
        Number of burn-in steps.
    sample_steps : int
        Number of sampling steps.
    thin : int
        Thinning interval; saves only every thin number of steps. Default
        is 1 (saves all steps, no thinning).
    pool : pool object
        Pool object for parallelization. Can be instance of emcee.utils.MPIPool
        or multiprocessing.pool (or other pool-like object implementing map
        method).
    betas : np.array
        Array containing the values to use for beta = 1/temperature.
    pos : numpy.array
        Matrix of initial positions for the chain. If None (default) random
        positions are chosen from the prior. Assigning a position allows
        previously run chains to be extended.
    random_state : random state for Mersenne Twister PRNG
        The random state to use to initialize the sampler's pseudo-random
        number generator. Can be used to continue runs from previous ones.
    """
    if pos_filename is None:
        warnings.warn('pos_filename was not specified, will not be able to '
                      'save intermediate burn-in positions.')

    # Initialize the parameter array with initial values (in log10 units)
    # Number of parameters to estimate
    ndim = (len(gf.builder.global_params) +
            (len(gf.data) * len(gf.builder.local_params)))
    # Initialize the walkers with starting positions drawn from the priors
    # Note that the priors are in log10 scale already, so they don't
    # need to be transformed here
    if pos is None:
        p0 = np.zeros((ntemps, nwalkers, ndim))
        for temp_ix in range(ntemps):
            for walk_ix in range(nwalkers):
                for p_ix in range(ndim):
                    p0[temp_ix, walk_ix, p_ix] = gf.priors[p_ix].random()
    else:
        p0 = pos

    # Create the sampler
    sampler = emcee.PTSampler(ntemps, nwalkers, ndim, likelihood, prior,
                              loglargs=[gf], logpargs=[gf], pool=pool,
                              betas=betas)

    # The PTSampler is implemented as a generator, so it is called in a for
    # loop
    # If we're not doing any burn-in, jump straight to sampling
    if burn_steps == 0:
        print "Main sampling..."
        nstep = 0
        for p, lnprob, lnlike in sampler.sample(p0, iterations=sample_steps,
                                            thin=thin):
            if nstep % 10 == 0:
                print "nstep %d of %d, MAP: %f" % (nstep, sample_steps,
                                                   np.max(lnprob[0]))
            nstep +=1
    # Otherwise, do the burn-in first
    else:
        print "Burn in sampling..."
        nstep = 0
        done = False
        last_ti = None
        print_interval = 1
        cur_start_position = p0
        abs_tol = 3.0 # The maximum allowable difference for convergence
        rel_tol = 0.1 # The fraction of the err allowable for convergence
        # Run the chain for rounds of convergence_interval steps; at the end
        # of each round, check for convergence. If converged, go on to main
        # sampling. If not, reinitialize sampler and run again. Running the
        # sampler in small rounds like this reduces the amount of memory
        # needed to just enough to store the chain for 1 round.
        while not done:
            if (burn_steps - nstep) < convergence_interval:
                num_iterations = burn_steps - nstep
            else:
                num_iterations = convergence_interval
            # Don't run again if we've already run more than the prescribed
            # number of burn in steps
            if num_iterations <= 0:
                break;
            for p, lnprob, lnlike in sampler.sample(cur_start_position,
                            iterations=num_iterations, storechain=True):
                # Increase the step counter by 1 since by the time we've gotten
                # here we've run an iteration of the sampler
                nstep += 1
                if nstep % print_interval == 0:
                    print("nstep %d of %d, MAP: %f, mean post %f" %
                         (nstep, burn_steps, np.max(lnprob[0]),
                          np.mean(lnprob[0])))
                    print sampler.tswap_acceptance_fraction
                    # Save the current position
                    if pos_filename is not None:
                        with open(pos_filename, 'w') as f:
                            rs = np.random.get_state()
                            cPickle.dump((p, rs), f)
                nstep += 1
            # If this is our first time checking convergence, set the TI
            # value and continue
            if last_ti is None:
                (last_ti, last_ti_err) = \
                            sampler.thermodynamic_integration_log_evidence()
                print "-- Initial TI value: %f, %f" % (last_ti, last_ti_err)
            else:
                (cur_ti, cur_ti_err) = \
                            sampler.thermodynamic_integration_log_evidence()
                diff = np.abs(last_ti - cur_ti)
                print("-- Last: %f, %f Current: %f Diff: %f" %
                      (last_ti, last_ti_err, cur_ti, diff))
                # Check for convergence
                if diff < abs_tol and cur_ti_err < abs_tol and \
                   last_ti_err < abs_tol and \
                   diff < (last_ti_err * rel_tol) and \
                   check_convergence_corr(sampler, 0, None,
                                          pval_threshold=0.001):
                    print "-- Converged!"
                    done = True
                else:
                    last_ti = cur_ti
                    last_ti_err = cur_ti_err
            # Reset the initial position to our last position
            cur_start_position = p
            # Reset the sampler
            sampler.reset()

        print "Main sampling..."
        nstep = 0
        for p, lnprob, lnlike in sampler.sample(p, lnprob0=lnprob,
                                     lnlike0=lnlike,
                                     iterations=sample_steps, thin=thin):
            # Increase the step counter by 1 since by the time we've gotten here
            # we've run an iteration of the sampler
            nstep += 1
            if nstep % 5 == 0:
                print "nstep %d of %d, MAP: %f, mean post %f" % \
                     (nstep, burn_steps, np.max(lnprob[0]), np.mean(lnprob[0]))
                print sampler.tswap_acceptance_fraction
            # Save the current position
            if pos_filename is not None:
                with open(pos_filename, 'w') as f:
                    rs = np.random.get_state()
                    cPickle.dump((p, rs), f)
        (final_ti, final_ti_err) = \
                sampler.thermodynamic_integration_log_evidence()
        print("-- Final TI: %f, %f --" % (final_ti, final_ti_err))
    # Close the pool!
    if pool is not None:
        pool.close()

    print "Done sampling."
    return sampler

def geweke_convergence(sampler, burn_frac=0.1, sample_frac=0.5,
                      p_threshold=0.05):
    # Define a few useful numbers
    ntemps = sampler.chain.shape[0]
    num_steps = sampler.chain.shape[2]
    num_params = sampler.chain.shape[3]
    burn_ubound = int(burn_frac * num_steps)
    nburn = burn_ubound - 1
    sample_lbound = int(sample_frac * num_steps)
    nsample = num_steps - sample_lbound
    converged = True # We've converged until we prove otherwise
    # Calculation of the test statistic
    def T_func(xm, xvar, nx, ym, yvar, ny):
        T = (xm - ym) / np.sqrt(xvar / float(nx) + yvar / float(ny))
        return np.abs(T)
    # Iterate over all temperatures
    for temp_ix in range(0, ntemps):
        print "Temp %d" % temp_ix
        # Iterate over all parameters
        for p_ix in range(num_params):
            # Get burn-in and sample steps
            burn_steps = sampler.chain[temp_ix, :, 0:burn_ubound, p_ix]
            sample_steps = sampler.chain[temp_ix, :, sample_lbound:, p_ix]
            # Calculate means and variances
            burn_mean = np.mean(burn_steps)
            sample_mean = np.mean(sample_steps)
            burn_var = np.var(burn_steps)
            sample_var = np.var(sample_steps)
            T = T_func(sample_mean, sample_var, nsample,
                       burn_mean, burn_var, nburn)
            if T < p_threshold:
                plt.ion()
                plt.figure()
                plt.plot(sampler.chain[temp_ix, :, :, p_ix].T, alpha=0.1)
                plt.plot(np.mean(sampler.chain[temp_ix, :, :, p_ix], axis=0))
                print("T = %f: not converged!" % T)
                converged = False
            else:
                print("T = %f" % T)
            # Now, check for convergence of product
            # g = (x_i - X)(y_i - Y)
            for q_ix in range(p_ix):
                print "param %d by %d" % (p_ix, q_ix)
                p_steps = sampler.chain[temp_ix, :, :, p_ix]
                q_steps = sampler.chain[temp_ix, :, :, q_ix]
                g_steps = (p_steps - np.mean(p_steps)) * \
                          (q_steps - np.mean(q_steps))
                g_burn_steps = g_steps[:, 0:burn_ubound]
                g_sample_steps = g_steps[:, sample_lbound:]
                g_burn_mean = np.mean(g_burn_steps)
                g_sample_mean = np.mean(g_sample_steps)
                g_burn_var = np.var(g_burn_steps)
                g_sample_var = np.var(g_sample_steps)
                T = T_func(g_sample_mean, g_sample_var, nsample,
                           g_burn_mean, g_burn_var, nburn)
                if T < p_threshold:
                    plt.ion()
                    plt.figure()
                    plt.plot(sampler.chain[temp_ix, :, :, p_ix].T, alpha=0.1)
                    plt.plot(np.mean(sampler.chain[temp_ix, :, :, p_ix],
                             axis=0))
                    print("T = %f: not converged!" % T)
                    converged = False
                else:
                    print("T = %f" % T)
    return converged

def check_convergence_corr(sampler, start_step, end_step, pval_threshold=0.2):
    # Only do the check if we've got enough steps
    if sampler.lnprobability is None:
        print "sampler.lnprobability is None."
        return

    lnpost = sampler.lnprobability[:, :, start_step:end_step]
    ntemps = lnpost.shape[0]
    nwalkers = lnpost.shape[1]
    nsteps = lnpost.shape[2]

    step_indices = np.repeat(np.arange(0, nsteps), nwalkers)
    pass_arr = np.zeros(ntemps)
    for temp_ix in range(0, ntemps):
        pts = lnpost[temp_ix].flatten(order='F')
        (rval, pval) = pearsonr(step_indices, pts)
        print "Temp %d: r = %f, p = %f" % (temp_ix, rval, pval)
        if rval < 0 or pval > pval_threshold:
            pass_arr[temp_ix] = True
    passes = np.all(pass_arr)
    print "Passes: ", passes
    return np.all(pass_arr)


def global_fit_from_args(args):
    ### DATA
    # Import the module containing the data
    data_args = args['data']
    __import__(data_args['module'])
    data_module = sys.modules[data_args['module']]
    # Get the relevant variables from the module containing the data
    data_var = data_module.__dict__[data_args['data_var']]
    data_sigma_var = data_module.__dict__[data_args['data_sigma_var']]
    time_var = data_module.__dict__[data_args['time_var']]
    # Get the name of the variable containing the initial conditions vector,
    # which may not exist
    ic_var_name = data_args['initial_condition_var']
    if ic_var_name is None:
        ic_var = None
    else:
        ic_var = data_module.__dict__[ic_var_name]

    ### MODEL
    # Call the appropriate model-building macro
    # If there is a multiconf attribute, that trumps any other attribute
    # and determines that this is a multiconf model
    if 'multiconf' in args['model']:
        bd = multiconf.Builder()
        # Number of confs and whether the model is reversible
        try:
            num_confs = int(args['model']['multiconf'])
            reversible = False
        # If it's not an int, assume it's two-element list with a flag
        # specifying a reversible model, indicated by 'rev'
        except TypeError:
            num_confs = int(args['model']['multiconf'][0])
            rev = args['model']['multiconf'][1]
            if rev == 'rev':
                reversible = True
            else:
                raise Exception('Unknown multiconf model flag %s' % rev)
        norm_data = args['model']['normalized_nbd_data']
        nbd_ubound = data_module.__dict__[data_args['nbd_ubound']]
        nbd_lbound = data_module.__dict__[data_args['nbd_lbound']]
        nbd_f0 = data_module.__dict__[data_args['nbd_f0']]
        bd.build_model_multiconf(num_confs, nbd_f0, nbd_lbound, nbd_ubound,
                                 normalized_data=norm_data,
                                 reversible=reversible)
    # Check the builder: one_cpt or lipo_sites
    elif 'builder' in args['model'] and \
         args['model']['builder'] == 'lipo_sites':
        bd = lipo_sites.Builder()
        bd.build_model_from_dict(args['model'])
    # If the builder is one_cpt or is not specified
    else:
        bd = one_cpt.Builder()
        bd.build_model_from_dict(args['model'])

    # Set the initial conditions
    for ic_name, ic_value in args['global_initial_conditions'].iteritems():
        bd.model.parameters[ic_name].value = ic_value

    ### PARAMETERS TO FIT
    if args['global_params'] == 'all':
        bd.global_params = bd.estimate_params
        bd.local_params = []
    else:
        bd.global_params = [bd.model.parameters[p_name]
                            for p_name in args['global_params']]
        bd.local_params = [bd.model.parameters[p_name]
                           for p_name in args['local_params']]

    local_ic_name = args['local_initial_condition']
    if ic_var is None or local_ic_name is None:
        params = None
    else:
        params = {local_ic_name: ic_var}

    # Create the global fit instance
    gf = GlobalFit(bd, time_var, data_var, data_sigma_var, params,
                   args['model_observable'])
    return gf
