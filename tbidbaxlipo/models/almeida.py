#from tbidbaxlipo.models import core, one_cpt
from pysb import builder
from bayessb.priors import Normal

class Builder(builder.Builder):

    def declare_monomers(self):
        """Declares signatures for tBid and Bax."""
        self.monomer('Bax',
                        ['loc', 'dye'],
                        {'loc': ['c', 'm', 'i', 'p'],
                         'dye': ['none', 'f', 'e']})
        self.monomer('Vesicles', ['dye'], {'dye': ['f', 'e']})

    def __init__(self, params_dict=None, nbd_sites=None, scaling_factor=None):
        # Sets self.model = Model(), and self.param_dict
        builder.Builder.__init__(self, params_dict=params_dict)

        self.declare_monomers()

        # INITIAL CONDITIONS
        self.parameter('Bax_0', 100)
        self.parameter('Vesicles_0', 0.7)

        Bax = self['Bax']
        Vesicles = self['Vesicles']

        self.initial(Bax(loc='c', dye='none'), self['Bax_0'])
        self.initial(Vesicles(dye='f'), self['Vesicles_0'])

        # OBSERVABLES
        self.observable('eVes', Vesicles(dye='e'))
        self.expression('Release', self['eVes'] / self['Vesicles_0'])

    def build_model_original(self):
        """The original model published for Cecropin A."""
        # Declare parameters
        kon = self.parameter('kon', 1e-4)
        koff = self.parameter('koff', 1e1)
        k1 = self.parameter('k1', 1e-2)
        k2 = self.parameter('k2', 1e1)
        keflx = self.parameter('keflx', 1e2)
        Bax = self['Bax']
        Vesicles = self['Vesicles']

        # Rules
        # Protein binds to membranes
        self.rule('Bax_binds_full',
                  Bax(loc='c', dye='none') + Vesicles(dye='f') >>
                  Bax(loc='m', dye='f') + Vesicles(dye='f'),
                  kon)
        self.rule('Bax_unbinds_full',
                  Bax(loc='m', dye='f') >> Bax(loc='c', dye='none'), koff)
        self.rule('Bax_binds_empty',
                  Bax(loc='c', dye='none') + Vesicles(dye='e') >>
                  Bax(loc='m', dye='e') + Vesicles(dye='e'),
                  kon)
        self.rule('Bax_unbinds_empty',
                  Bax(loc='m', dye='e') >> Bax(loc='c', dye='none'), koff)

        # On full membranes, Bax enters pore state
        self.rule('Bax_forms_pore',
                  Bax(loc='m', dye='f') >> Bax(loc='p', dye='f'), k1)
        # Bax goes from pore state to empty membrane bound state
        self.rule('Bax_full_to_empty',
                  Bax(loc='p', dye='f') >> Bax(loc='m', dye='e'), k2)
        # Dye release proportional to P*
        self.rule('Dye_release',
                  Vesicles(dye='f') + Bax(loc='p', dye='f') >>
                  Vesicles(dye='e') + Bax(loc='p', dye='f'),
                  keflx)


    def build_model_bax_bound(self):
        """The original model published for Cecropin A."""
        # Declare parameters
        kon = self.parameter('kon', 1e-3)
        koff = self.parameter('koff', 1e1)
        kp_fwd = self.parameter('kp_fwd', 1e-2)
        kp_rev = self.parameter('kp_rev', 0)
        k_pore = self.parameter('k_pore', 1e1)
        keflx = self.parameter('keflx', 1e2)
        Bax = self['Bax']
        Vesicles = self['Vesicles']

        # Rules
        # Protein binds to membranes
        self.rule('Bax_binds_full',
                  Bax(loc='c', dye='none') + Vesicles(dye='f') >>
                  Bax(loc='m', dye='f') + Vesicles(dye='f'),
                  kon)
        self.rule('Bax_unbinds_full',
                  Bax(loc='m', dye='f') >> Bax(loc='c', dye='none'), koff)
        self.rule('Bax_binds_empty',
                  Bax(loc='c', dye='none') + Vesicles(dye='e') >>
                  Bax(loc='m', dye='e') + Vesicles(dye='e'),
                  kon)
        self.rule('Bax_unbinds_empty',
                  Bax(loc='m', dye='e') >> Bax(loc='c', dye='none'), koff)

        # On membranes, Bax enters pore state
        self.rule('Bax_forms_pore_fwd',
                  Bax(loc='m') >> Bax(loc='p'), kp_fwd)
        # Bax goes from pore state to empty membrane bound state
        self.rule('Bax_full_to_empty',
                  Bax(loc='p', dye='f') >> Bax(loc='p', dye='e'), k_pore)
        # Bax goes from pore state back to m state and can unbind
        self.rule('Bax_forms_pore_rev',
                  Bax(loc='p', dye='e') >> Bax(loc='m', dye='e'), kp_rev)
        # Dye release proportional to P*
        self.rule('Dye_release',
                  Vesicles(dye='f') + Bax(loc='p', dye='f') >>
                  Vesicles(dye='e') + Bax(loc='p', dye='f'),
                  keflx)

        self.parameter('pore_scaling', 0.7)
        self.expression('Release', self['eVes'] / self['pore_scaling'])

    def build_model_bax_two_state(self):
        """The original model published for Cecropin A."""
        # Declare parameters
        kon = self.parameter('kon', 1e-3)
        koff = self.parameter('koff', 1e-1)
        k_act_fwd = self.parameter('k_act_fwd', 1.8e0)
        kp_fwd = self.parameter('kp_fwd', 5e-3)
        kp_rev = self.parameter('kp_rev', 0)
        k_pore = self.parameter('k_pore', 3e3)
        keflx = self.parameter('keflx', 1e2)
        Bax = self['Bax']
        Vesicles = self['Vesicles']

        # Rules
        # Protein binds to membranes
        self.rule('Bax_binds_full',
                  Bax(loc='c', dye='none') + Vesicles(dye='f') >>
                  Bax(loc='m', dye='f') + Vesicles(dye='f'),
                  kon)
        self.rule('Bax_unbinds_full',
                  Bax(loc='m', dye='f') >> Bax(loc='c', dye='none'), koff)
        self.rule('Bax_binds_empty',
                  Bax(loc='c', dye='none') + Vesicles(dye='e') >>
                  Bax(loc='m', dye='e') + Vesicles(dye='e'),
                  kon)
        self.rule('Bax_unbinds_empty',
                  Bax(loc='m', dye='e') >> Bax(loc='c', dye='none'), koff)

        # On membranes, Bax enters activated state
        self.rule('Bax_insertion_fwd',
                  Bax(loc='m') >> Bax(loc='i'), k_act_fwd)
        # Activatd state goes to pore state
        self.rule('Bax_forms_pore_fwd',
                  Bax(loc='i') >> Bax(loc='p'), kp_fwd)
        # Bax goes from pore state to empty membrane bound state
        self.rule('Bax_full_to_empty',
                  Bax(loc='p', dye='f') >> Bax(loc='p', dye='e'), k_pore)
        # Bax goes from pore state back to m state and can unbind
        self.rule('Bax_forms_pore_rev',
                  Bax(loc='p', dye='e') >> Bax(loc='m', dye='e'), kp_rev)
        # Dye release proportional to P*
        self.rule('Dye_release',
                  Vesicles(dye='f') + Bax(loc='p', dye='f') >>
                  Vesicles(dye='e') + Bax(loc='p', dye='f'),
                  keflx)

        self.observable('iBax_obs', self['Bax'](loc='i'))
        self.expression('iBax_frac', self['iBax_obs'] / self['Bax_0'])
        self.observable('cBax_obs', self['Bax'](loc='c'))
        self.expression('cBax_frac', self['cBax_obs'] / self['Bax_0'])


if __name__ == '__main__':
    from pysb.bng import generate_equations
    from pysb.integrate import Solver
    from tbidbaxlipo.util import fitting
    from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
    import numpy as np
    from matplotlib import pyplot as plt
    import sys

    plt.ion()

    builder = Builder()
    #builder.build_model()
    #builder.build_model_bax_bound()
    builder.build_model_bax_two_state()

    activator = 'Bid'
    nbd_site = 'WT'
    rep_index = 1
    rt = df[(activator, 'Release', nbd_site, rep_index, 'TIME')].values
    ry = df[(activator, 'Release', nbd_site, rep_index, 'VALUE')].values / 100.

    plt.figure()
    rp = -np.log(1 - ry)
    plt.plot(rt, rp)
    plt.title('Pores')

    rpd = np.diff(rp)
    rtd = rt[1:]
    plt.figure()
    plt.plot(rtd, rpd)
    plt.title('Pore derivative')
    sys.exit()

    a = fitting.Parameter(1e-3)
    b = fitting.Parameter(1e-4)
    c = fitting.Parameter(1e-4)
    def gompertz(t):
        return a() * np.exp(-b() * np.exp(-c()* t))
    fitting.fit(gompertz, [a, b, c], rp, rt)
    plt.plot(rt, gompertz(rt))

    plt.figure()
    plt.plot(rt, ry, color='r')
    plt.title('Release')




    #t = np.linspace(0, 10000, 10000)
    t = rt
    solver = Solver(builder.model, t)
    solver.run()
    plt.figure()
    plt.plot(rt, ry, color='r')
    plt.plot(t, solver.yexpr['Release'])
    plt.plot(t, solver.yexpr['iBax_frac'])
    plt.plot(t, solver.yexpr['cBax_frac'])

    pysb_fit = fitting.fit_pysb_builder(builder, 'Release', rt, ry)
    plt.figure()
    plt.plot(rt, ry, color='r')
    plt.plot(rt, pysb_fit.ypred)

