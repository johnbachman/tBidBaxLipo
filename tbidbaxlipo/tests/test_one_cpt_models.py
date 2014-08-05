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
    plt.figure()
    plt.plot(t, x)
