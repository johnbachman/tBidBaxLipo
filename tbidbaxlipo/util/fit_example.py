import fitting
import numpy as np

t = np.linspace(0, 100, 50)
data = 3*t + 4

m = fitting.Parameter(2)
b = fitting.Parameter(5)

def fit_func(x): return (m()*x + b())

fitting.fit(fit_func, [m, b], data, t)

print "m: %f" % m()
print "b: %f" % b()


