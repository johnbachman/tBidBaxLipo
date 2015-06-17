import yaml
import sys
from tbidbaxlipo.util import emcee_fit
import numpy as np
from scipy import optimize
from pyDOE import lhs # Latin hypercube sampling
import cPickle

def deterministic_fit(gf, p0, posterior, method='Nelder-Mead', bounds=None):
    res = None
    if method == 'Nelder-Mead':
        # Downhill Simplex or "Nelder-Mead" algorithm
        print("Method: Nelder-Mead/Simplex (default)")
        res = optimize.fmin(posterior, p0, args=(gf,), maxiter=10000,
                            maxfun=10000, full_output=True)
    elif method == 'BFGS':
        print("Method: BFGS")
        res = optimize.fmin_l_bfgs_b(posterior, p0, args=(gf,), fprime=None,
                                     approx_grad=True, bounds=bounds)
    else:
        raise ValueError("Unknown method %s" % method)

    return res

def generate_latin_hypercube(gf, num_samples, basename):
    # The number of parameters
    ndim = (len(gf.builder.global_params) +
            (len(gf.data) * len(gf.builder.local_params)))
    # Generate a latin hypercube of the parameters
    lh = lhs(ndim, samples=num_samples, criterion='center')
    # For each sample...
    for samp_ix in range(num_samples):
        # ...initialize the vector of initial values
        p0 = np.zeros(ndim)
        # For each parameter...
        for p_ix in range(ndim):
            # ...get the prior...
            pr = gf.priors[p_ix]
            # ...and use the inverse CDF of the prior distribution to convert
            # the percentile value on [0, 1] to an initial value
            percentile = lh[samp_ix, p_ix]
            p0[p_ix] = pr.inverse_cdf(percentile)
        print("Saving position: %s" % p0)
        filename = '%s.%d.lhs' % (basename, samp_ix)
        with open(filename, 'w') as f:
            cPickle.dump(p0, f)
    return

def fit(gf, p0):
    # The size of our initial position vector should always match the number
    # of parameters we're fitting!
    assert len(gf.priors) == p0.shape[0]
    # Specify the amount by which we tweak the lower and upper bounds to prevent
    # log probabilities of -inf (for use by L-BFGS-B algorithm)
    epsilon = 0.001

    # Before fitting, get the list of bounds for constrained optimization
    bounds = []
    # For each parameter...
    for p_ix in range(p0.shape[0]):
        # ...get the prior...
        pr = gf.priors[p_ix]
        # If the prior has hard bounds, add to the list of parameter bounds
        if hasattr(pr, 'lower_bound') and hasattr(pr, 'upper_bound'):
            bounds_tup = (pr.lower_bound + epsilon, pr.upper_bound - epsilon)
            bounds.append(bounds_tup)
        else:
            bounds.append((None, None))
            print bounds

    # Run the fit!
    res = deterministic_fit(gf, p0, emcee_fit.negative_posterior,
                            method='Nelder-Mead', bounds=bounds)

if __name__ == '__main__':
    import os.path

    # Check arguments
    usage_msg = "Usage:\n"
    usage_msg += "%s hypercube num_samples yaml_file\n" % \
                  os.path.basename(__file__)
    usage_msg += "%s fit yaml_file lhs_file\n" % \
                  os.path.basename(__file__)

    if len(sys.argv) < 4:
        print usage_msg
        sys.exit(1)
    # TODO consolidate args
    elif sys.argv[1] == 'hypercube' and len(sys.argv) != 4:
        print usage_msg
        sys.exit(1)
    elif sys.argv[1] == 'fit' and len(sys.argv) != 4:
        print usage_msg
        sys.exit(1)

    # Seed the random number generator to create the Latin Hypercube
    np.random.seed(1)

    # Get the command type

    # HYPERCUBE
    if sys.argv[1] == 'hypercube':
        # ...load the args for the data/model from the YAML file
        with open(sys.argv[3]) as yaml_file:
            print("Loading YAML file %s" % sys.argv[3])
            args = yaml.load(yaml_file)
        # Get the GlobalFit object from the args
        gf = emcee_fit.global_fit_from_args(args)
        # Generate the hypercube
        num_samples = int(sys.argv[2])
        basename = sys.argv[3]
        generate_latin_hypercube(gf, num_samples, basename)
    # FIT
    elif sys.argv[1] == 'fit':
        # ...load the args for the data/model from the YAML file
        with open(sys.argv[2]) as yaml_file:
            print("Loading YAML file %s" % sys.argv[2])
            args = yaml.load(yaml_file)
        # Get the GlobalFit object from the args
        gf = emcee_fit.global_fit_from_args(args)
        # Load the initial position
        p0_filename = sys.argv[3]
        with open(p0_filename) as p0_file:
            print("Loading initial position file %s" % sys.argv[2])
            p0 = cPickle.load(p0_file)
        # Run the fit
        res = fit(gf, p0)
        # Save the results
        result_filename = '%s.detfit' % os.path.splitext(p0_filename)[0]
        with open(result_filename, 'w') as result_file:
            cPickle.dump((gf, res), result_file)
    # (error)
    else:
        print(usage_msg)
        sys.exit(1)

    # Pick initial guess from the priors
    #for p_ix in range(ndim):
    #    p0[p_ix] = gf.priors[p_ix].random()

