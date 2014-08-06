from tbidbaxlipo.models.one_cpt import Builder
from pysb.integrate import odesolve
import numpy as np
from matplotlib import pyplot as plt

t = np.linspace(0, 100, 100)
plt.ion()

def test_build_model_t():
    b = Builder()
    b.build_model_t()
    x = odesolve(b.model, t)
    plot_model(t, x)

def test_build_model_ta():
    b = Builder()
    b.build_model_ta()
    x = odesolve(b.model, t)
    plot_model(t, x)

def plot_model(t, x):
    plt.figure()
    for name in x.dtype.names:
        if not name.startswith('__'):
            plt.plot(t, x[name], label=name)
    plt.legend(loc='upper right')
