import numpy as np
from tbidbaxlipo.models import core
from tbidbaxlipo.priors import Uniform, UniformLinear

class Builder(core.Builder):

    def __init__(self, params_dict=None):
        core.Builder.__init__(self, params_dict=params_dict)

    def build_model_multiconf(self, num_confs, c0_scaling, nbd_lbound=None,
                              nbd_ubound=None, normalized_data=False,
                              reversible=False):
        if num_confs < 2:
            raise ValueError('There must be a minimum of two conformations.')

        self.num_confs = num_confs
        self.reversible = reversible

        # Initialize monomer and initial condition
        Bax = self.monomer('Bax', ['conf'],
                           {'conf': ['c%d' % i for i in range(num_confs)]})
        Bax_0 = self.parameter('Bax_0', 1, prior=None)
        self.initial(Bax(conf='c0'), Bax_0)

        # Scaling for initial conformation
        scaling = self.parameter('c0_scaling', c0_scaling, prior=None)
        obs = self.observable('Bax_c0', Bax(conf='c0'))
        sympy_expr = (scaling * obs)

        # Set the bounds for the scaling parameters
        if normalized_data:
            scaling_prior = UniformLinear(-1, 1)
        else:
            if nbd_lbound is None or nbd_ubound is None:
                raise ValueError("If NBD data is not normalized, upper and "
                                 "lower bounds for the scaling parameters "
                                 "must be explicitly specified.")
            scaling_prior = UniformLinear(np.log10(nbd_lbound),
                                          np.log10(nbd_ubound))

        # Rules for transitions between other conformations
        for i in range(num_confs-1):
            rate = self.parameter('c%d_to_c%d_k' % (i, i+1), 1e-3,
                                  prior=Uniform(-6, -1))
            scaling = self.parameter('c%d_scaling' % (i+1), 1,
                                     prior=scaling_prior)

            self.rule('c%d_to_c%d' % (i, i+1),
                      Bax(conf='c%d' % i) >> Bax(conf='c%d' % (i+1)), rate)
            if reversible:
                rate = self.parameter('c%d_to_c%d_k' % (i+1, i), 1e-3,
                                      prior=Uniform(-6, -1))
                self.rule('c%d_to_c%d' % (i+1, i),
                          Bax(conf='c%d' % (i+1)) >> Bax(conf='c%d' % i), rate)

            obs = self.observable('Bax_c%d' % (i+1), Bax(conf='c%d' % (i+1)))

            sympy_expr += (scaling * obs)

        # The expression mapping to our experimental observable
        self.expression('NBD', sympy_expr)

        # Set the model name
        self.model.name = "%dconfs" % num_confs

        # Set the model nbd func
        if self.num_confs == 2 or (self.num_confs == 3 and reversible == False):
            self.set_nbd_func()


    def set_nbd_func(self):
        """Assigns a function to self.formula that, when called after setting
        the parameter values of self.model, returns the timecourse for the
        given parameters."""
        if self.num_confs == 2 and self.reversible == False:
            def nbd_func(t):
                Bax_0 = self['Bax_0'].value
                c0 = Bax_0 * np.exp(-self['c0_to_c1_k'].value * t)
                nbd = (self['c0_scaling'].value * c0 +
                       self['c1_scaling'].value * (Bax_0 - c0))
                return nbd
        elif self.num_confs == 2 and self.reversible == True:
            def nbd_func(t):
                kf = self['c0_to_c1_k'].value
                kr = self['c1_to_c0_k'].value
                Bax_0 = self['Bax_0'].value
                c0eq = (kr * Bax_0) / float(kf + kr)
                c1eq = (kf * Bax_0) / float(kf + kr)
                c0 = c1eq * np.exp(-(kf+kr)*t) + c0eq
                nbd = (self['c0_scaling'].value * c0 +
                       self['c1_scaling'].value * (Bax_0 - c0))
                return nbd
        elif self.num_confs == 3 and self.reversible == False:
            k1 = self['c0_to_c1_k'].value
            k2 = self['c1_to_c2_k'].value
            Bax_0 = self['Bax_0'].value
            # First, the case where k1 and k2 are equal (need to treat
            # separately to avoid divide by 0 errors)
            if k1 == k2:
                def nbd_func(t):
                    c0 = Bax_0 * np.exp(-k1 * t)
                    c1 = Bax_0 * t * k1 * np.exp(-k1 * t)
                    nbd = (self['c0_scaling'].value * c0 +
                           self['c1_scaling'].value * c1 +
                           self['c2_scaling'].value * (Bax_0 - c0 - c1))
                    return nbd
            # The typical case, where k1 and k2 are not equal
            else:
                def nbd_func(t):
                    Bax_0 = self['Bax_0'].value
                    c0 = Bax_0 * np.exp(-k1 * t)
                    c1 = (((np.exp(-k1*t) - np.exp(-k2*t)) * k1 * Bax_0) /
                          (k2 - k1))
                    nbd = (self['c0_scaling'].value * c0 +
                           self['c1_scaling'].value * c1 +
                           self['c2_scaling'].value * (Bax_0 - c0 - c1))
                    return nbd
        self.nbd_func = nbd_func





