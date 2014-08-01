from pysb import *
import numpy as np
from tbidbaxlipo.models import core
from bayessb.priors import Uniform, Normal
from pysb.integrate import Solver
from matplotlib import pyplot as plt

class Builder(core.Builder):

    def __init__(self, params_dict=None, nbd_sites=None):
        core.Builder.__init__(self, params_dict=params_dict)

        Plate = self.monomer('Plate', ['bax'])
        Bax = self.monomer('Bax',
                           ['plate', 'bleached', 'ins', 'ves', 'mem', 'nbd'],
                           {'bleached': ['n','y'], 'ins': ['n', 'y'], 'nbd': ['y', 'n'],
                            'mem': ['n', 'y']})
        Vesicles = self.monomer('Vesicles', ['bax'])

        # Initial conditions for Bax
        Bax_NBD_0 = self.parameter('Bax_NBD_0', 184, prior=None)
        Bax_unlab_0 = self.parameter('Bax_unlab_0', 0, prior=None)

        # Initial conditions
        self.initial(Bax(bleached='n', plate=None, ves=None, mem='n', ins='n', nbd='y'),
                     Bax_NBD_0)
        self.initial(Bax(bleached='n', plate=None, ves=None, mem='n', ins='n', nbd='n'),
                     Bax_unlab_0)

    def background_processes(self):
        Bax = self['Bax']
        Plate = self['Plate']

        Plate_0 = self.parameter('Plate_0', 20, prior=Uniform(0, 2))
        k_bleach = self.parameter('k_bleach', 1e-5, prior=Normal(-5, 2))
        k_plate_bind = self.parameter('k_plate_bind', 1e-5, prior=Normal(-5, 2))
        cBax_NBD = self.parameter('cBax_NBD', 5.5/185.,
                                prior=Normal(np.log10(5.5/185), np.log10((0.15/185.)**2)))

        # Rules
        self.rule('Bax_sticks_to_plate',
             Bax(plate=None, ves=None, ins='n') + Plate(bax=None) >>
             Bax(plate=1, ves=None, ins='n') % Plate(bax=1),
             k_plate_bind)
        self.rule('Bax_photobleaches',
             Bax(bleached='n') >> Bax(bleached='y'),
             k_bleach)

        self.initial(Plate(bax=None), Plate_0)

    def competitive_lipo_binding(self):
        Bax = self['Bax']
        Vesicles = self['Vesicles']
        cBax_NBD = self['cBax_NBD']

        # Parameters for insertion mechanism
        iBax_NBD = self.parameter('iBax_NBD', 4,
                                  prior=Uniform(np.log10(0.9), np.log10(6)))
        kf_Bax_binding = self.parameter('kf_Bax_binding', 1e-4, prior=Normal(-4, 3))
        kr_Bax_binding = self.parameter('kr_Bax_binding', 1e-1, prior=Normal(-1, 4))
        k_Bax_insertion = self.parameter('k_Bax_insertion', 1e-4, prior=Normal(-4, 3))
        Vesicles_0 = self.parameter('Vesicles_0', 0.015, prior=None)

        self.rule('Bax_binds_vesicles',
                  Bax(plate=None, ves=None, ins='n') + Vesicles(bax=None) <>
                  Bax(plate=None, ves=1, ins='n') % Vesicles(bax=1),
                  kf_Bax_binding, kr_Bax_binding)
        self.rule('Bax_inserts_into_vesicles',
                  Bax(plate=None, ves=1, ins='n') % Vesicles(bax=1) >>
                  Bax(plate=None, ves=None, ins='y') + Vesicles(bax=None),
                  k_Bax_insertion)

        self.initial(Vesicles(bax=None), self['Vesicles_0'])

        # Outputs: observable and expression for NBD fluorescence
        cBax_ = self.observable('cBax_', Bax(bleached='n', plate=None, ins='n', nbd='y'))
        iBax_ = self.observable('iBax_', Bax(bleached='n', plate=None, ins='y', nbd='y'))
        self.expression('NBD', cBax_ * cBax_NBD + (iBax_ * iBax_NBD * cBax_NBD))

    def noncompetitive_lipo_binding(self):
        Bax = self['Bax']
        Vesicles = self['Vesicles']
        cBax_NBD = self['cBax_NBD']

        # Parameters for insertion mechanism
        iBax_NBD = self.parameter('iBax_NBD', 4,
                                  prior=Uniform(np.log10(0.9), np.log10(6)))
        kf_Bax_binding = self.parameter('kf_Bax_binding', 1e-4, prior=Normal(-4, 3))
        kr_Bax_binding = self.parameter('kr_Bax_binding', 1e-1, prior=Normal(-1, 4))
        k_Bax_insertion = self.parameter('k_Bax_insertion', 1e-4, prior=Normal(-4, 3))
        Vesicles_0 = self.parameter('Vesicles_0', 0.015, prior=None)

        self.rule('Bax_binds_vesicles',
                  Bax(plate=None, mem='n', ins='n') + Vesicles(bax=None) >>
                  Bax(plate=None, mem='y', ins='n') + Vesicles(bax=None),
                  kf_Bax_binding)
        self.rule('Bax_unbinds_vesicles',
                  Bax(plate=None, mem='y', ins='n') >>
                  Bax(plate=None, mem='n', ins='n'),
                  kr_Bax_binding)

        self.rule('Bax_inserts_into_vesicles',
                  Bax(plate=None, mem='y', ins='n') >>
                  Bax(plate=None, mem='y', ins='y'),
                  k_Bax_insertion)

        self.initial(Vesicles(bax=None), self['Vesicles_0'])

        # Outputs: observable and expression for NBD fluorescence
        cBax_ = self.observable('cBax_', Bax(bleached='n', plate=None, ins='n', nbd='y'))
        iBax_ = self.observable('iBax_', Bax(bleached='n', plate=None, ins='y', nbd='y'))
        self.expression('NBD', cBax_ * cBax_NBD + (iBax_ * iBax_NBD * cBax_NBD))

    def build_competitive_model(self):
        self.background_processes()
        self.competitive_lipo_binding()

    def build_noncompetitive_model(self):
        self.background_processes()
        self.noncompetitive_lipo_binding()

    def bax_titration_experiment(self, t, ydata):
        self['Bax_NBD_0'].value = 200.
        self['Vesicles_0'].value = 1.9
        bax_concs = np.array([0, 10, 20, 40, 80, 160, 320, 640])

        plt.figure()
        for bax_conc in bax_concs:
            self['Bax_unlab_0'].value = bax_conc
            baseline = self['cBax_NBD'].value * self['Bax_NBD_0'].value
            solver = Solver(self.model, t)
            solver.run()
            plt.plot(t, solver.yexpr['NBD'] / baseline, label='Bax %d' % bax_conc)
            print "baseline = %f" % baseline
        plt.legend(loc='upper left')

        plt.plot(t, ydata / (self['cBax_NBD'].value * 185))

if __name__ == '__main__':
    plt.ion()

    from tbidbaxlipo.plots.layout_140318 import bgsub_wells, TIME, VALUE
    bg = bgsub_wells['A4']
    t = bg[TIME]
    ydata = bg[VALUE]

    max_posterior_params = {
        'iBax_NBD': 4.3240275339272909,
        'kf_Bax_binding': 0.006306510357230467,
        'kr_Bax_binding': 0.33527514631180588,
        'k_Bax_insertion': 0.017148023284087228,
        'Plate_0': 17.509258995395331,
        'cBax_NBD': ydata[0]/185.,
        'k_bleach': 1.0442068851723714e-05,
        'k_plate_bind': 1.212354253716292e-05,
    }
    b = Builder(params_dict=max_posterior_params)
    b.build_competitive_model()
    b.bax_titration_experiment(t, ydata)
