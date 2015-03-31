from tbidbaxlipo.util import error_propagation as ep
from nose.tools import ok_
import numpy as np

def test_calc_ratio_mean_sd_num():
    numer_mean = 1
    denom_mean = 4
    numer_sd =0.25
    denom_sd = 0.8
    result = ep.calc_ratio_mean_sd(numer_mean, numer_sd, denom_mean, denom_sd)
    # Check that return value is an iterable (should be a tuple) with two
    # entries, each a number
    ok_(isinstance(result, tuple))
    ok_(len(result) == 2)

def test_calc_ratio_mean_sd_array():
    numer_mean = np.array([1, 2, 3])
    denom_mean = np.array([2, 4, 6])
    numer_sd = np.array([0.25, 0.3, 0.35])
    denom_sd = np.array([0.5, 0.6, 0.7])
    result = ep.calc_ratio_mean_sd(numer_mean, numer_sd, denom_mean, denom_sd)
    # Check that return value is an iterable (should be a tuple) with two
    # entries, each a number
    ok_(isinstance(result, tuple))
    (ratio_mean, ratio_sd) = result
    ok_(isinstance(ratio_mean, np.ndarray))
    ok_(isinstance(ratio_sd, np.ndarray))
    ok_(len(ratio_mean.shape) == 1)
    ok_(len(ratio_sd.shape) == 1)
    ok_(len(ratio_mean) == len(ratio_sd))
    ok_(len(ratio_mean) == len(numer_mean))
