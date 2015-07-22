"""
Implementation of prior distributions for MCMC fitting.

These classes provide the ability to calculate the probability of given values
according the appropriate PDF, and also allow for sampling of random numbers
from these distributions (which is useful for choosing randomized initial
values for an MCMC walk, for example).

I chose to write these simple classes rather than use the distributions in
scipy.stats because according to a few empirical tests using the %time command
in ipython, I found that making calculations from a uniform PDF (which would be
required at every MCMC iteration) was approximately two orders of magnitude
faster using this approach rather than the classes in scipy.stats.
"""

import numpy as np
import scipy.stats

class Uniform():
    """A uniform prior distribution.

    Parameters
    ----------
    lower_bound : number
        The lower bound of the interval.
    upper_bound : number
        The upper bound of the interval.
    """
    def __init__(self, lower_bound, upper_bound):
        self.lower_bound = float(lower_bound)
        self.upper_bound = float(upper_bound)
        self.p = 1. / (upper_bound - lower_bound)

    def pdf(self, num):
        """Returns the negative log probability of the given number.

        The negative log probability is given by:

        infinity, if num is outside of the bounds of the interval;
        log(1 / (upper - lower)), otherwise.
        """
        if num < self.lower_bound or num > self.upper_bound:
            return np.inf
        else:
            return -np.log(self.p)

    def random(self):
        """Get a random sample from the uniform distribution."""
        return np.random.uniform(low=self.lower_bound, high=self.upper_bound)

    def inverse_cdf(self, percentile):
        """Get the value associated with the percentile probability.

        Also known as the percent point function. A wrapper around
        scipy.stats.uniform.ppf().
        """
        scale = self.upper_bound - self.lower_bound
        return scipy.stats.uniform.ppf(percentile, loc=self.lower_bound,
                                       scale=scale)

class UniformLinear(Uniform):
    """A uniform prior distribution that inverts a log10-transformation.

    Parameters
    ----------
    lower_bound : number
        The lower bound of the interval.
    upper_bound : number
        The upper bound of the interval.
    """
    def __init__(self, lower_bound, upper_bound):
        self.linear_lower_bound = 10 ** float(lower_bound)
        self.linear_upper_bound = 10 ** float(upper_bound)
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.p = 1. / (self.linear_upper_bound - self.linear_lower_bound)

    def pdf(self, num):
        """Returns the negative log probability of the given number.

        The negative log probability is given by:

        infinity, if num is outside of the bounds of the interval;
        -log(1 / (upper - lower)), otherwise.
        """
        num = 10 ** num
        if num < self.linear_lower_bound or num > self.linear_upper_bound:
            return np.inf
        else:
            return -np.log(self.p)

    def random(self):
        """Get a random sample from the (non-log) uniform distribution."""
        rand_sample = np.random.uniform(low=self.linear_lower_bound,
                                        high=self.linear_upper_bound)
        return np.log10(rand_sample)

    def inverse_cdf(self, percentile):
        """Get the value associated with the percentile probability.

        Also known as the percent point function. A wrapper around
        scipy.stats.uniform.ppf().
        """
        scale = self.linear_upper_bound - self.linear_lower_bound
        value = scipy.stats.uniform.ppf(percentile, loc=self.linear_lower_bound,
                                        scale=scale)
        return np.log10(value)

class Normal():
    def __init__(self, mean, stdev):
        self.mean = float(mean)
        self.stdev = float(stdev)

    def pdf(self, num):
        """Returns the negative log probability of the given number.

        The negative log probability is given by:

        (num - mean)^2 / (2 * stdev ** 2)
        """
        return ((num - self.mean)**2) / (2. * self.stdev ** 2)

    def random(self):
        """Get a random sample from the normal distribution."""
        return (self.stdev * np.random.randn()) + self.mean

    def inverse_cdf(self, percentile):
        """Get the value associated with the percentile probability.

        Also known as the percent point function. A wrapper around
        scipy.stats.norm.ppf().
        """
        return scipy.stats.norm.ppf(percentile, loc=self.mean,
                                    scale=self.stdev)

