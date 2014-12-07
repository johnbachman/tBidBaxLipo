from pysb import *
import pysb.builder
from bayessb.priors import Normal, Uniform
import numpy as np
from pysb.integrate import Solver
from matplotlib import pyplot as plt

class Builder(pysb.builder.Builder):
    def __init__(self, params_dict=None):
        # Sets self.model = Model(), and self.param_dict
        super(Builder, self).__init__(params_dict=params_dict)

        Bid = self.monomer('Bid', ['dye', 'b'], {'dye': ['y', 'n']})
        Lipos = self.monomer('Lipos', ['b'])

        Lipos_0 = self.parameter('Lipos_0', 10., prior=Uniform(0, 4)) # nM
        kf = self.parameter('kf', 1e-3, prior=Uniform(-6, -3))
        kr = self.parameter('kr', 1e-3, prior=Uniform(-4, 1))
        # Initial conditions, don't estimate
        Bid_lab_0 = self.parameter('Bid_lab_0', 10., prior=None)
        Bid_unlab_0 = self.parameter('Bid_unlab_0', 0., prior=None)

        FRET_efficiency = self.parameter('FRET_efficiency', 35, prior=Uniform(0,2))

        # Set initial conditions
        self.initial(Bid(dye='y', b=None), Bid_lab_0)
        self.initial(Bid(dye='n', b=None), Bid_unlab_0)
        self.initial(Lipos(b=None), Lipos_0)
        # Binding rule
        self.rule('Bid_binds_Lipos',
                  Bid(b=None) + Lipos(b=None) <> Bid(b=1) % Lipos(b=1),
                  kf, kr)
        Bid_lab_bound = self.observable('Bid_lab_bound',
                                         Bid(b=1, dye='y') % Lipos(b=1))
        self.expression('FRET', (Bid_lab_bound / Bid_lab_0) * FRET_efficiency)

if __name__ == '__main__':
    bd = Builder()
    t = np.linspace(0, 8000, 1000)
    sol = Solver(bd.model, t)
    sol.run()
    plt.ion()
    plt.figure()
    plt.plot(t, sol.yexpr['FRET'])
