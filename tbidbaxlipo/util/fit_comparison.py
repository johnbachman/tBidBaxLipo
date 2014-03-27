import numpy as np
import pymc
from tbidbaxlipo.util import fitting
from matplotlib import pyplot as plt

"""
# Bayesian parameter estimation for biochemical kinetics

Application I: '''Estimating uncertainty in parameter values.'''
Motivation: Is liposome background release rate invariant with concentration?
Motivation: What is variation expected for background release?

Gradient-based
--------------
In the gradient-based method, you fit all the replicates at the same time and
the parameters minimize the error across the replicates.  End up with a maximum
likelihood (true in the nonlinear case? I think so, if the errors are presumed
in N(0,1)) estimate of the parameter values.  Moreover, the covariances of the
parameters can be determined from the covariance matrix, and the residuals can
be used to estimate the confidence intervals on the parameters. So you do get
some estimate of the uncertainty in your fits.

***TODO*** I wonder, would scaling the objective function by the SDs of the
points have the same effects as fitting all the replicates simultaneously?
Should try this for the fits of a single point.

***TODO*** So need modified fit function to fit replicates--that takes a matrix
of y values and a single vector of input values, and turns it into a matched
set of vectors to do the subtraction.

Bayesian
--------

In the Bayesian method, you could either fit the average, with the associated
SDs, or you could fit the whole set of replicates simultaneously as done for
the gradient-based method. Again, is this the same, when normal errors are
assumed?

Running MCMC on an average of even 5 or 10 reps produces an estimate of the
uncertainty in the parameters. Naturally, the uncertainty goes down as the
number of replicates goes up (and hence the standard error goes down). This
allows for good estimates of the means underlying the various distributions.

One advantage is that the MCMC gives you the full posterior distribution which
is useful if the posterior distributions are not themselves normal. You can
then look at the order statistics of the distribution, which is useful if the
posterior distribution has significant skew. Should see if this is the case in
the liposome release case.

Hierarchical Fitting
====================

So, is this adequate? Is it important to know the variance in the replicate
distribution or is it adequate to simply know the uncertainty associated with
the average?  I suppose the issue is that if you don't know the spread
associated with single replicates, you don't know the range of possibilities
for the background associated with a particular well--you will simply subtract
use the mean value, with any associated uncertainty.  It doesn't make sense
that the variation in background release that may be affecting your
experimental condition should be affected by the number of replicates you do to
estimate the mean. In fact, for this reason it would seem that what you really
want to do is estimate the mean and variance associated with the individual
wells.

Could 
"""


# Samples
num_pts = 10000
sigma = 5
pts = sigma * np.random.randn(num_pts)

mu = fitting.Parameter(np.random.randn(1))

def mu_fit(t):
    return mu()

result = fitting.fit(mu_fit, [mu], pts, 1)

error_variance = np.var(pts - mu_fit(1))

plt.ion()
plt.figure()
plt.hist(pts, bins=10)
plt.vlines(mu(), 0, plt.gca().get_ylim()[1], linewidth=2, color='r')

std_err = np.std(pts) / np.sqrt(num_pts)

cov_x = result[1]
print "SD from cov_x: %f" % np.sqrt(cov_x[0] * error_variance)
print "SD from sample: %f" % std_err
print "True stderr: %f" % (sigma/np.sqrt(num_pts))

# So, to get the true covariance of the parameters, I have to fit

