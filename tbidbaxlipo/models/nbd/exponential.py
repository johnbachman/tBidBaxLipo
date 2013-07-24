from tbidbaxlipo.models import core
from bayessb.priors import Uniform

class Builder(core.Builder):

    def __init__(self, params_dict=None):
        core.Builder.__init__(self, params_dict=params_dict)

    def build_model_exponential(self, num_exponentials, f0):
        if num_exponentials < 1:
            raise ValueError('There must be a minimum of one exponential.')

        self.num_exponentials = num_exponentials

        # All initial conditions are 1
        Bax_0 = self.parameter('Bax_0', 1, estimate=False)

        # Scaling for initial signal
        f0_param = self.parameter('f0', f0, estimate=False)
        sympy_expr = f0_param

        # Initialize monomer and initial condition
        for i in range(num_exponentials):
            m = self.monomer('Bax_exp%d' % i, ['conf'], {'conf': ['c0', 'c1']})
            self.initial(m(conf='c0'), Bax_0)
            rate = self.parameter('Bax_exp%d_k' % i, 1e-3,
                                  prior=Uniform(-6, 1))
            scaling = self.parameter('Bax_exp%d_fmax' % i, 1e3,
                                     prior=Uniform(-1, 5))
            self.rule('Bax_exp%d_rule' % i, m(conf='c0') >> m(conf='c1'), rate)
            obs = self.observable('Bax_exp%d_obs' % i, m(conf='c1'))
            sympy_expr += (obs*scaling)

        # The expression mapping to our experimental observable
        self.expression('NBD', sympy_expr)

        # Set the name of the model
        self.model.name = "%dexp" % num_exponentials
