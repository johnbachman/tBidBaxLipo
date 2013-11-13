import numpy as np

def calc_ratio_sd(numer_mean, numer_sd, denom_mean, denom_sd):
    """Calculates the variance of a ratio of two normal distributions with
    the given means and standard deviations."""
    numer_samples = numer_mean + (numer_sd * np.random.randn(1000))
    denom_samples = denom_mean + (denom_sd * np.random.randn(1000))
    ratio_samples = numer_samples / denom_samples
    return np.std(ratio_samples)

