from tbidbaxlipo.models.one_cpt import Builder
from pysb.integrate import odesolve
import numpy as np
from matplotlib import pyplot as plt

t = np.linspace(0, 8000, 100)
plt.ion()

def test_build_models():
    model_names = ['t', 'ta', 'tai', 'taid', 'taidt', 'tair', 'taird',
                   'tairdt', 'tad', 'tadt', 'tar', 'tard', 'tardt']
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
