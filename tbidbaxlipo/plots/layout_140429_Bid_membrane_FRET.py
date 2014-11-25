import numpy as np
from matplotlib import pyplot as plt
import tbidbaxlipo.data
import sys
import os
from tbidbaxlipo.util import fitting

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
data_file = os.path.abspath(os.path.join(data_path,
                                        '140429_Bid568_DiDlipo_FRET.csv'))
fret_data = np.loadtxt(data_file, delimiter=',')

plt.ion()
plt.close('all')
#plt.figure()
#plt.plot(r, cf(r), marker='o')

bid_concs = fret_data[:,0] + 10.
fret = fret_data[:,1]

def fit_basic_binding_curve():
    kd = fitting.Parameter(6.) # kd for binding lipo sites
    f = fitting.Parameter(30.) # FRET efficiency
    f0 = fitting.Parameter(10.) # baseline quenching

    def fit_func(bid_conc):
        frac_bound = bid_conc / (bid_conc + kd())
        frac_dye = 10. / bid_conc
        frac_dye_bound = frac_bound * frac_dye
        return f() * frac_dye_bound + f0()

    fitting.fit(fit_func, [kd, f], fret, bid_concs)

    plt.figure('Bid FRET')
    plt.plot(np.log10(bid_concs), fret, marker='o')
    plt.plot(np.log10(bid_concs), fit_func(bid_concs))
    print "Kd: %f" % kd()
    print "f: %f" % f()
    print "f0: %f" % f0()

def fit_stoichiometric_comp():
    kd = fitting.Parameter(6.) # kd for binding lipo sites
    f = fitting.Parameter(30.) # FRET efficiency
    f0 = fitting.Parameter(10.) # baseline quenching

    def fit_func(bid_conc):
        frac_bound = bid_conc / (bid_conc + kd())
        frac_dye = 10. / bid_conc
        frac_dye_bound = frac_bound * frac_dye
        return f() * frac_dye_bound + f0()

    fitting.fit(fit_func, [kd, f], fret, bid_concs)

    plt.figure('Bid FRET')
    plt.plot(np.log10(bid_concs), fret, marker='o')
    plt.plot(np.log10(bid_concs), fit_func(bid_concs))
    print "Kd: %f" % kd()
    print "f: %f" % f()
    print "f0: %f" % f0()




def fit_gouy_chapman():
    lipid_conc = 129.8 # uM; 0.1 (mg/mL); ~1.5 nM liposomes

    v = fitting.Parameter(6.) # charges per peptide chain
    b = fitting.Parameter(1.) # dimensionless
    pc = fitting.Parameter(1.) # partition coefficient
    f = fitting.Parameter(35.) # FRET efficiency
    f0 = fitting.Parameter(9.5) # baseline quenching
    r = np.logspace(-3, 0, 1000)

    def cf(r):
        return (r * alpha(r)) / pc()
    def alpha(r):
        log_alpha = 2 * v() * np.arcsinh(v() * b() * r)
        return np.exp(log_alpha)

    def fit_func(bid_concs):
        cb = r * lipid_conc
        cfr = cf(r)
        ctot = cb + cfr
        fracs_bound = cb / ctot
        interp_fracs = np.interp(bid_concs, ctot, fracs_bound)
        return f() * interp_fracs + f0()

    plt.figure()
    plt.plot(np.log10(bid_concs), fret, marker='o') # Data
    plt.plot(np.log10(bid_concs), fit_func(bid_concs))

    fitting.fit(fit_func, [v, b, pc, f, f0], fret, bid_concs)
    plt.figure()
    plt.plot(np.log10(bid_concs), fret, marker='o') # Data
    plt.plot(np.log10(bid_concs), fit_func(bid_concs))

    print f0()
    print f()
    print v()
    print b()
    print pc()
    import ipdb; ipdb.set_trace()


if __name__ == '__main__':
    fit_basic_binding_curve()
    fit_gouy_chapman()
    fit_no_competition()

