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

        plt.figure()
        plt.plot(t, self.obs_func(x))
        plt.show()

    def build_2_conf_model(self):
        self.build_model_multiconf(2)
        def obs_func(sim_data):
            return sim_data['Bax_c1']*self.model.parameters['c1_scaling'].value
        self.obs_func = obs_func
        self.model.name = '2_conf_model'

    def build_3_conf_model(self):
        self.build_model_multiconf(3)
        def obs_func(sim_data):
            return ((sim_data['Bax_c1'] *
                     self.model.parameters['c1_scaling'].value) +
                    (sim_data['Bax_c2'] *
                     self.model.parameters['c2_scaling'].value))
        self.obs_func = obs_func
        self.model.name = '3_conf_model'

    def build_4_conf_model(self):
        self.build_model_multiconf(4)
        def obs_func(sim_data):
            return ((sim_data['Bax_c1'] *
                     self.model.parameters['c1_scaling'].value) +
                    (sim_data['Bax_c2'] *
                     self.model.parameters['c2_scaling'].value) +
                    (sim_data['Bax_c3'] *
                     self.model.parameters['c3_scaling'].value))
        self.obs_func = obs_func
        self.model.name = '4_conf_model'


