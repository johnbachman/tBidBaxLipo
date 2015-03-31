from tbidbaxlipo.models import core
from bayessb.priors import Uniform

class Builder(core.Builder):

    def __init__(self, params_dict=None):
        core.Builder.__init__(self, params_dict=params_dict)

    def build_model_multiconf(self, num_confs, c0_scaling,
                              normalized_data=False, reversible=False):
        if num_confs < 2:
            raise ValueError('There must be a minimum of two conformations.')

        self.num_confs = num_confs

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
            scaling_prior = Uniform(-2, 2)
        else:
            scaling_prior = Uniform(1, 5)

        # Rules for transitions between other conformations
        for i in range(num_confs-1):
            rate = self.parameter('c%d_to_c%d_k' % (i, i+1), 1e-3,
                                  prior=Uniform(-6, -1))
            scaling = self.parameter('c%d_scaling' % (i+1), 1,
                                     prior=scaling_prior)

            self.rule('c%d_to_c%d' % (i, i+1),
                      Bax(conf='c%d' % i) >> Bax(conf='c%d' % (i+1)), rate)
            if reversible:
                rate = self.parameter('c%d_to_c%d_kr' % (i+1, i), 1e-3,
                                      prior=Uniform(-6, -1))
                self.rule('c%d_to_c%d' % (i+1, i),
                          Bax(conf='c%d' % (i+1)) >> Bax(conf='c%d' % i), rate)

            obs = self.observable('Bax_c%d' % (i+1), Bax(conf='c%d' % (i+1)))

            sympy_expr += (scaling * obs)

        # The expression mapping to our experimental observable
        self.expression('NBD', sympy_expr)

        # Set the model name
        self.model.name = "%dconfs" % num_confs
