from tbidbaxlipo.models import core, one_cpt
from bayessb.priors import Normal

class Builder(one_cpt.Builder):

    def declare_monomers(self):
        """Declares signatures for tBid and Bax."""
        self.monomer('tBid', ['bh3', 'loc'],
                {'loc': ['c', 'm']})
        self.monomer('Bax',
                        ['bh3', 'a6', 'loc', 'lipo'],
                        {'loc':  ['c', 'm', 'i', 'a', 'p', 'bh3expos']})
        self.monomer('Vesicles', ['bax', 'dye'],
                     {'dye': ['f', 'e']})

    def __init__(self, params_dict=None, nbd_sites=None, scaling_factor=None):
        # Sets self.model = Model(), and self.param_dict
        core.Builder.__init__(self, params_dict=params_dict)

        self.declare_monomers()

        # COMPARTMENTS
        solution = self.compartment('solution', dimension=3, parent=None)
        self.compartment('ves', dimension=2, parent=solution)

        # INITIAL CONDITIONS
        self.parameter('Vesicles_0', 1.9, prior=None)
        #self.parameter('tBid_0', 20, prior=None)
        self.parameter('Bax_0', 100, prior=None)

        tBid = self['tBid']
        Bax = self['Bax']
        Vesicles = self['Vesicles']

        #self.initial(tBid(loc='c', bh3=None) ** solution, self['tBid_0'])
        self.initial(Bax(loc='c', bh3=None, a6=None, lipo=None) ** solution,
                     self['Bax_0'])
        self.initial(Vesicles(bax=None, dye='f') ** solution, self['Vesicles_0'])

        # OBSERVABLES
        #self.observable('ctBid', tBid(loc='c'))
        #self.observable('mtBid', tBid(loc='m'))
        self.observable('cBax', Bax(loc='c'))
        self.observable('mBax', Bax(loc='m'))
        self.observable('iBax', Bax(loc='i'))
        self.observable('aBax', Bax(loc='a'))
        self.observable('pBax', Bax(loc='p'))
        self.observable('eVes', Vesicles(dye='e'))
        # self.observable('tBidBax', tBid(bh3=1) % Bax(bh3=1))
        # self.observable('Bax2', Bax(bh3=1) % Bax(bh3=1))
        #self.observable('Baxbh3', Bax(bh3=1))
        #self.observable('Bax4',
        #     MatchOnce(Bax(loc='i', bh3=1, a6=3) % Bax(loc='i', bh3=1, a6=4) % 
        #               Bax(loc='i', bh3=2, a6=3) % Bax(loc='i', bh3=2, a6=4)))
        # Pore formation
        #self.observable('pores', Pores())

        # SCALING PARAMETERS
        #if nbd_sites is not None:
        #    self.declare_nbd_scaling_parameters(nbd_sites)

    def build_model_nbd_terbium(self):
        # Bax gets to the membrane
        self.translocate_Bax()
        # Bax is activated
        self.basal_Bax_activation(reversible=True)
        # Get monomers from the model
        Bax = self['Bax']
        Vesicles = self['Vesicles']
        # Bax catalytically permeabilizes vesicles
        self.parameter('iBax_enz_release_k', 1e-4,
                       prior=Normal(-4, 3))
        self.rule('iBax_enz_release',
                  Bax(loc='i') + Vesicles(dye='f') >>
                  Bax(loc='i') + Vesicles(dye='e'),
                  self['iBax_enz_release_k'])
        # Bax converts to a pore-incompetent state
        # Like c1->c2
        self.parameter('iBax_to_pBax_kf', 1e-3,
                       prior=Normal(-3, 3))
        self.parameter('iBax_to_pBax_kr', 1e-4,
                       prior=Normal(-4, 3))
        self.rule('iBax_to_pBax',
                  Bax(loc='i') <> Bax(loc='p'),
                  self['iBax_to_pBax_kf'], self['iBax_to_pBax_kr'])
        # The observable: fraction released
        c0_scaling = self.parameter('c0_scaling', 1, prior=None)
        c1_scaling = self.parameter('c1_scaling', 1.5,
                                prior=Normal(0, 0.5))
        c2_scaling = self.parameter('c2_scaling', 2,
                                prior=Normal(0, 0.5))
        # Get observables
        self.expression('NBD', (c0_scaling * (self['cBax'] + self['mBax']) +
                               c1_scaling * self['iBax'] +
                               c2_scaling * self['pBax']) / self['Bax_0'])
        self.expression('Tb', self['eVes'] / self['Vesicles_0'])
        self.model.name = 'nbd_terbium'


if __name__ == '__main__':
    from pysb.integrate import Solver
    from matplotlib import pyplot as plt
    import numpy as np

    b = Builder()
    b.build_model_nbd_terbium()
    t = np.linspace(0, 1e5)
    s = Solver(b.model, t)
    s.run()
    plt.ion()
    plt.figure()
    plt.plot(t, s.yexpr['NBD'])
    plt.figure()
    plt.plot(t, s.yexpr['Tb'])
    
