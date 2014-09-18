import pysb.builder
from tbidbaxlipo.util import fitting
from tbidbaxlipo.data.parse_bid_bim_fret_nbd_release import df
from bayessb.priors import Normal
import numpy as np
from pysb.integrate import Solver
import numpy as np
from matplotlib import pyplot as plt

class Builder(pysb.builder.Builder):
    """
    Class/constructor documentation here.
    """
    def __init__(self, params_dict=None, nbd_sites=None):
        # Sets self.model = Model(), and self.param_dict
        super(Builder, self).__init__(params_dict=params_dict)
        self.cpt_list = ['ves']
        self.declare_components()

    def declare_components(self):
        """Declares signatures for tBid, Bax, and Vesicles.

        Must be called after the child class constructor has initialized the
        list of compartments.
        """
        tBid = self.monomer('tBid', ['bh3', 'conf', 'cpt'],
                     {'conf': ['aq', 'mem'],
                      'cpt':  ['sol'] + self.cpt_list})
        Bax = self.monomer('Bax', ['bh3', 'a6', 'lipo', 'conf', 'cpt', 'pore'],
                     {'conf':  ['aq', 'mem', 'ins'],
                      'cpt': ['sol'] + self.cpt_list,
                      'pore': ['y', 'n']})
        Vesicles = self.monomer('Vesicles', ['bax'])
        Pores = self.monomer('Pores', ['cpt'], {'cpt':self.cpt_list})
        self.parameter('tBid_0', 20)
        self.parameter('Bax_0', 100)

        self.initial(tBid(cpt='ves', conf='mem', bh3=None), self['tBid_0'])
        self.initial(Bax(cpt='ves', conf='mem', bh3=None, a6=None,
                         lipo=None, pore='n'),
                     self['Bax_0'])
        #self.initial(Vesicles(bax=None), self['Vesicles_0'])

        # OBSERVABLES
        self.observable('ctBid', tBid(conf='aq'))
        self.observable('mtBid', tBid(conf='mem'))

        self.observable('cBax', Bax(cpt='sol', conf='aq'))
        self.observable('mBax', Bax(conf='mem'))
        self.observable('iBax', Bax(conf='ins'))
        self.observable('mBax_free', Bax(conf='mem', bh3=None))
        self.observable('iBax_free', Bax(conf='ins', bh3=None))
        self.observable('tBidiBax', tBid(bh3=1) % Bax(bh3=1, conf='ins'))
        self.observable('tBidmBax', tBid(bh3=1) % Bax(bh3=1, conf='mem'))
        self.observable('Bax2', Bax(bh3=1) % Bax(bh3=1))
        self.observable('Bax4',
                        Bax(conf='ins', bh3=1, a6=3) %
                        Bax(conf='ins', bh3=1, a6=4) %
                        Bax(conf='ins', bh3=2, a6=3) %
                        Bax(conf='ins', bh3=2, a6=4))
        # Pore formation
        self.observable('pBax', Bax(pore='y'))
        self.observable('pores', Pores())

        for cpt_name in self.cpt_list:
            self.observable('mtBid_%s' % cpt_name,
                            tBid(conf='mem', cpt=cpt_name))
            self.observable('mBax_%s' % cpt_name,
                            Bax(conf='mem', cpt=cpt_name))
            self.observable('iBax_%s' % cpt_name,
                            Bax(conf='ins', cpt=cpt_name))
            self.observable('tBidBax_%s' % cpt_name,
                            tBid(bh3=1, cpt=cpt_name) %
                            Bax(bh3=1, cpt=cpt_name))
            self.observable('Bax2_%s' % cpt_name,
                            Bax(bh3=1, cpt=cpt_name) %
                            Bax(bh3=1, cpt=cpt_name))
            self.observable('Bax4_%s' % cpt_name,
                            Bax(conf='ins', bh3=1, a6=3, cpt=cpt_name) %
                            Bax(conf='ins', bh3=1, a6=4, cpt=cpt_name) %
                            Bax(conf='ins', bh3=2, a6=3, cpt=cpt_name) %
                            Bax(conf='ins', bh3=2, a6=4, cpt=cpt_name))
            # Pore formation
            self.observable('pBax_%s' % cpt_name, Bax(pore='y', cpt=cpt_name))
            self.observable('pores_%s' % cpt_name, Pores(cpt=cpt_name))

    def random_initial_values(self):
        initial_values_log = np.empty(len(self.priors))
        for i, prior in enumerate(self.priors):
            initial_values_log[i] = prior.random()
        return 10 ** initial_values_log

    def prior(self, mcmc, position):
        """Calculate the prior probability of the set of parameters.

        The builder maintains a list of prior objects that corresponds in its
        order to the list of parameters in the parameter list (it includes only
        the parameters to be estimated).  To calculate the prior, the method
        iterates over the objects in self.priors and for each one gets the
        negative log probability of that position in parameter space. The
        negative log probabilities are then summed (equivalent to the product
        of the individual priors) and then returned by the function.
        """
        prior_prob = 0
        for i, prior in enumerate(self.priors):
            prior_prob += prior.pdf(position[i])
        return prior_prob

    def build_model_fret1(self):
        tBid = self['tBid']
        Bax = self['Bax']

        p1 = self.parameter('tBid_mBax_kf', 0.01, prior=Normal(-2, 2))
        p2 = self.parameter('tBid_mBax_kr', 0.01, prior=Normal(-2, 2))
        p3 = self.parameter('tBidmBax_tBidiBax_kf', 0.001, prior=Normal(-3, 2))
        p4 = self.parameter('tBidmBax_tBidiBax_kr', 0.0001, prior=Normal(-4, 2))
        p5 = self.parameter('tBidiBax_tBid_iBax_kf', 0.01, prior=Normal(-2, 2))
        p6 = self.parameter('tBidiBax_tBid_Bax_kr', 0.001, prior=Normal(-3, 2))

        self.rule('tBid_binds_mBax',
                  tBid(bh3=None) + Bax(bh3=None, conf='mem') <>
                  tBid(bh3=1) % Bax(bh3=1, conf='mem'),
                  p1, p2)

        self.rule('tBidmBax_to_tBidiBax',
                  tBid(bh3=1) % Bax(bh3=1, conf='mem') <>
                  tBid(bh3=1) % Bax(bh3=1, conf='ins'),
                  p3, p4)

        self.rule('tBidiBax_to_tBid_iBax',
                  tBid(bh3=1) % Bax(bh3=1, conf='mem') <>
                  tBid(bh3=None) + Bax(bh3=None, conf='ins'),
                  p5, p6)

        c1_nbd = self.parameter('c1_nbd', 2, prior=Normal(0, 1))
        c2_nbd = self.parameter('c2_nbd', 3, prior=Normal(0, 1))
        c3_nbd = self.parameter('c3_nbd', 5., prior=Normal(0, 1))
        c1_fret = self.parameter('c1_fret', 30, prior=Normal(1, 1))
        c2_fret = self.parameter('c2_fret', 20, prior=Normal(1, 1))

        self.expression('FRET', (self['tBidmBax'] * c1_fret +
                                 self['tBidiBax'] * c2_fret)
                                / float(self['tBid_0'].value))
        self.expression('NBD', (self['mBax_free'] +
                                self['tBidmBax'] * c1_nbd +
                                self['tBidiBax'] * c2_nbd +
                                self['iBax_free'] * c3_nbd)
                               / float(self['Bax_0'].value))

if __name__ == '__main__':
    plt.ion()

    bd = Builder()
    bd.build_model_fret1()
    t = np.linspace(0, 4000)
    s = Solver(bd.model, t)
    s.run()
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.plot(t, s.yexpr['NBD'])
    plt.title('NBD')
    plt.subplot(1, 2, 2)
    plt.plot(t, s.yexpr['FRET'])
    plt.title('FRET')

    t = df[('Bid', 'NBD', '126', 1, 'TIME')]
    ny = df[('Bid', 'NBD', '126', 1, 'VALUE')]
    fy = df[('Bid', 'FRET', '126', 1, 'VALUE')]
    plt.figure()
    plt.plot(t, ny)
    plt.plot(t, fy)
