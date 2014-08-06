from tbidbaxlipo.models.one_cpt import Builder
from pysb.integrate import odesolve
import numpy as np
from matplotlib import pyplot as plt

t = np.linspace(0, 10000, 100)
plt.ion()

def test_build_models():
    model_names = [#'t', 'ta', 'tai', 'taid', 'taidt', 'tair', 'taird',
                   #'tairdt', 'tad', 'tadt', 'tar', 'tard', 'tardt',
                    #---
                   #'bax_heat',
                   #'bax_heat_aggregation',
                   #'bax_heat_reversible',
                    #---
                   #'bax_heat_reversible_aggregation',
                   #'bax_heat_dimer',
                   #'bax_heat_dimer_reversible',
                   #'bax_heat_auto1',
                   #'bax_heat_auto2',
                   #'bax_heat_bh3_auto2',
                   #'bax_heat_auto1_reversible_activation',
                   #'bax_heat_auto2_reversible_activation',
                   #'bax_heat_auto1_reversible',
                   #'bax_heat_auto2_reversible',
                   #'bax_heat_auto1_dimer',
                   #'bax_heat_auto2_dimer',
                   #'bax_heat_auto1_dimer_reversible',
                   #'bax_heat_auto2_dimer_reversible',
                   #'bax_heat_bh3_exposure_auto2',
                   #'bax_schwarz',
                   #'bax_schwarz_reversible',
                   #'bax_schwarz_dimer',
                   #'bax_schwarz_dimer_reversible',
                   #'bax_schwarz_tetramer_reversible',
    ]

    for model_name in model_names:
        print "Running model %s" % model_name
        b = Builder()
        eval("b.build_model_%s()" % model_name)
        x = odesolve(b.model, t)
        plot_model(t, x, 'build_model_%s' % model_name)

def plot_model(t, x, title):
    plt.figure()
    for name in x.dtype.names:
        if not name.startswith('__'):
            plt.plot(t, x[name], label=name)
    plt.title(title)
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.legend(loc='upper right')
