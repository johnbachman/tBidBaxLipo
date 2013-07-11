
from pysb import *
from pysb.integrate import Solver, odesolve
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.models import core

class Builder(core.Builder):

    def __init__(self, params_dict=None):
        core.Builder.__init__(self, params_dict=params_dict)

    def build_model_multiconf(self, num_confs):
        if num_confs < 2:
            raise ValueError('There must be a minimum of two conformations.')

        self.num_confs = num_confs

        # Initialize monomers and observable
        Bax = self.monomer('Bax', ['conf'],
                           {'conf': ['c%d' % i for i in range(num_confs)]})
        Bax_0 = self.parameter('Bax_0', 1)
        self.initial(Bax(conf='c0'), Bax_0)
        self.observable('Bax_c0', Bax(conf='c0'))

        for i in range(num_confs-1):
            rate = self.parameter('c%d_to_c%d_k' % (i, i+1), 1e-3)
            scaling = self.parameter('c%d_scaling' % (i+1), 1)

            self.rule('c%d_to_c%d' % (i, i+1),
                      Bax(conf='c%d' % i) >> Bax(conf='c%d' % (i+1)), rate)
            self.observable('Bax_c%d' % (i+1), Bax(conf='c%d' % (i+1)))
        # observable should be a product of c0 and c1 and scaling factors

    def run(self, tmax=10):
        t = np.linspace(0, tmax, 1000)
        x = odesolve(self.model, t)
        plt.ion()
        plt.figure()
        for i in range(self.num_confs):
            plt.plot(t, x['Bax_c%d' % i], label='Bax_c%d' % i)
        plt.legend(loc='upper right')
        plt.show()

        #plt.figure()
        #obs = x['Bax_c0']*1 + x['Bax_c1']*0 + x['Bax_c2']*1 + x['Bax_c3']*0
        #plt.plot(t, obs)
        #plt.show()


