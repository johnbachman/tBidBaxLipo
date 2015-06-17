import yaml
import sys
from tbidbaxlipo.util import emcee_fit
import numpy as np
from scipy import optimize
from pyDOE import lhs # Latin hypercube sampling
import cPickle
import glob

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
        filename = '%s.%dof%d.lhs' % (basename, samp_ix+1, num_samples)
        with open(filename, 'w') as f:
            cPickle.dump(p0, f)
    return

def fit(gf, p0_filename):
    # Load the initial position from the pickle file
    with open(p0_filename) as p0_file:
        print("Loading initial position file %s" % p0_filename)
        p0 = cPickle.load(p0_file)
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
    # Save the results
    result_filename = '%s.detfit' % os.path.splitext(p0_filename)[0]
    with open(result_filename, 'w') as result_file:
        cPickle.dump((gf, res), result_file)
    # Return the results for interactive
    return res

def load_globalfit_from_file(yaml_filename):
    with open(yaml_filename) as yaml_file:
        print("Loading YAML file %s" % yaml_filename)
        args = yaml.load(yaml_file)
    # Get the GlobalFit object from the args
    gf = emcee_fit.global_fit_from_args(args)
    return gf

if __name__ == '__main__':
    import os.path
    this_file = os.path.basename(__file__)

    # Check arguments
    usage_msg = "Usage:\n"
    usage_msg += "%s hypercube num_samples yaml_file\n" % this_file
    usage_msg += "%s fit yaml_file lhs_file\n" % this_file
    usage_msg += "%s submit [none|bsub|qsub] yaml_file\n" % this_file

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

    # HYPERCUBE ------------------------------------------------
    if sys.argv[1] == 'hypercube':
        # Load the args for the data/model from the YAML file
        gf = load_globalfit_from_file(sys.argv[3])
        # Generate the hypercube
        num_samples = int(sys.argv[2])
        basename = sys.argv[3]
        generate_latin_hypercube(gf, num_samples, basename)
    # FIT ------------------------------------------------------
    elif sys.argv[1] == 'fit':
        # Load the args for the data/model from the YAML file
        gf = load_globalfit_from_file(sys.argv[2])
        # Load the initial position
        p0_filename = sys.argv[3]
        # Run the fit
        res = fit(gf, p0_filename)
    # SUBMIT ---------------------------------------------------
    # det_fit.py submit none yaml_filename
    elif sys.argv[1] == 'submit':
        scheduler = sys.argv[2]
        # No job scheduler--just run each of the files in a loop
        if scheduler == 'none':
            # Load the args for the data/model from the YAML file
            yaml_filename = sys.argv[3]
            gf = load_globalfit_from_file(yaml_filename)
            # Get the list of position files
            p0_files = glob.glob(r'%s.*of*.lhs' % yaml_filename)
            if not p0_files:
                raise Exception("No p0 files found!")
            # Iterate over the position files, loading and running each
            for p0_filename in p0_files:
                fit(gf, p0_filename)
        elif scheduler == 'qsub':
            # Load the args for the data/model from the YAML file
            yaml_filename = sys.argv[3]
            gf = load_globalfit_from_file(yaml_filename)
            # Get the list of position files
            p0_files = glob.glob(r'%s.*of*.lhs' % yaml_filename)
            if not p0_files:
                raise Exception("No p0 files found!")
            # Iterate over the position files, loading and running each
            for p0_filename in p0_files:
                p0_basename = os.path.splitext(p0_filename)[0]
                err_filename = '%s.err' % p0_basename
                out_filename = '%s.out' % p0_basename
                qsub_args = ['qsub', '-b', 'y', '-cwd', '-V', '-o',
                             out_filename, '-e', err_filename,
                             'python', '-m', 'tbidbaxlipo.pt.det_fit', 'fit',
                             yaml_filename, p0_filename]
                print ' '.join(qsub_args)
        else:
            raise ValueError()
    else:
        print(usage_msg)
        sys.exit(1)

    # Pick initial guess from the priors
    #for p_ix in range(ndim):
    #    p0[p_ix] = gf.priors[p_ix].random()

