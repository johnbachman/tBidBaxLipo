from pysb import *
from tbidbaxlipo.util import fitting
import numpy as np
from matplotlib import pyplot as plt
from pysb.integrate import Solver
import pysb.builder

plt.ion()

class Builder(pysb.builder.Builder):
    def __init__(self, params_dict=None):
        super(Builder, self).__init__(params_dict=params_dict)
        self.build_model()

    def build_model(self):

        A = self.monomer('A', ['b'])
        B = self.monomer('B', ['b'])

        kf = self.parameter('kf', 1e-3)
        kr = self.parameter('kr', 1e-3)
        A_0 = self.parameter('A_0', 100)
        B_0 = self.parameter('B_0', 100)

        self.initial(A(b=None), A_0)
        self.initial(B(b=None), B_0)

        self.rule('A_binds_B', A(b=None) + B(b=None) <> A(b=1) % B(b=1), kf, kr)

        self.observable('AB', A(b=1) % B(b=1))

        self.global_params = [kr]
        self.local_params = [kf]

# Generate some synthetic data
num_timecourses = 10
data = []
t = np.linspace(0, 60, 100)
B_initial_vect = np.logspace(0, 3, num_timecourses)
seed = 1

plt.figure()

bd = Builder()

for i in range(num_timecourses):
    s = Solver(bd.model, t)
    bd.model.parameters['B_0'].value = B_initial_vect[i]
    s.run()
    ysim = s.yobs['AB']
    plt.plot(t, ysim, color='b')

    # Add error to the underlying data
    sigma = 0.01
    ydata = ysim * (np.random.randn(len(ysim)) * sigma + 1)
    plt.plot(t, ydata, color='r')

    data.append(ydata)

params = {'B_0': B_initial_vect}

gf = fitting.GlobalFit(bd, t, data, params, 'AB', obs_type='Observable')
result = gf.fit(method='Nelder-Mead')

gf.plot_func(result.x)

