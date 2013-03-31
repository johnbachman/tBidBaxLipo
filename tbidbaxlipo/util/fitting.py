from numpy import shape, arange, sum
from scipy import optimize

### STUFF FOR PARAMETER FITTING
# from http://www.scipy.org/Cookbook/FittingData#head-5eba0779a34c07f5a596bbcf99dbc7886eac18e5
class Parameter:
    def __init__(self, value):
            self.value = value
    def set(self, value):
            self.value = value
    def __call__(self):
            return self.value

def fit(function, parameters, y, x = None, maxfev=100000):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        err = y - function(x)
        return err

    if x is None: x = arange(y.shape[0])
    p = [param() for param in parameters]
    result = optimize.leastsq(f, p, ftol=1e-12, xtol=1e-12, maxfev=maxfev)

def mse(function, y, x):
    err = y - function(x)
    n = err.shape[0]
    sq_err = err**2
    return sum(sq_err) / n

def residuals(function, y, x):
    err = y - function(x)
    return err

