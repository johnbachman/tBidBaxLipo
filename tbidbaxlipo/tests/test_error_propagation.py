from tbidbaxlipo.util import error_propagation as ep
from nose.tools import ok_

def test_calc_ratio_mean_sd_num():
    numer_mean = 1
    denom_mean = 4
    numer_sd =0.25
    denom_sd = 0.8
    result = ep.calc_ratio_mean_sd(numer_mean, numer_sd, denom_mean, denom_sd)
    # Check that return value is an iterable (should be a tuple) with two
    # entries, each a number
    ok_(type(result) is tuple)
    ok_(len(result) == 2)
