import numpy as np

def calc_ratio_sd(numer_mean, numer_sd, denom_mean, denom_sd, num_samples=1000):
    """Calculates the variance of a ratio of two normal distributions with
    the given means and standard deviations."""
    numer_samples = numer_mean + (numer_sd * np.random.randn(num_samples))
    denom_samples = denom_mean + (denom_sd * np.random.randn(num_samples))
    ratio_samples = numer_samples / denom_samples
    return np.std(ratio_samples)

def calc_ratio_mean_sd(numer_mean, numer_sd, denom_mean, denom_sd,
                       num_samples=10000):
    """Calculates the variance of a ratio of two normal distributions with
    the given means and standard deviations."""
    numer_samples = numer_mean + (numer_sd * np.random.randn(num_samples))
    denom_samples = denom_mean + (denom_sd * np.random.randn(num_samples))
    ratio_samples = numer_samples / denom_samples
    return (np.mean(ratio_samples), np.std(ratio_samples))
