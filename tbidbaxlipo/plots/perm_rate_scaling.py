"""Analysis of how the empirically observed rate of permeabiliation (as
measured by exponential k/Fmax plots) scales with both one-compartment and
multi-compartment models of permeabilization."""

from tbidbaxlipo.models import one_cpt
import numpy as np
from matplotlib import pyplot as plt
from pysb.integrate import Solver
from tbidbaxlipo.plots.titration_fits import OneExpFmax, TwoExpNoT

plt.ion()

bax_concs = np.logspace(-1, 2, 30)
all_params = []

for bax_conc in bax_concs:
    #params_dict = {'Bax_0':bax_conc, 'Vesicles_0':1., 'Bax_transloc_kr':0.}
    params_dict = {'Bax_0':bax_conc, 'Vesicles_0':1.}
    bd = one_cpt.Builder(params_dict=params_dict)
    bd.build_model_bax_schwarz_irreversible()

    t = np.linspace(0, 10000)
    s = Solver(bd.model, t)
    s.run()

    Ves_0 = float(bd.model.parameters['Vesicles_0'].value)
    dr = 1 - np.exp(-(s.yobs['pores']/Ves_0))

    plt.figure('Timecourses')
    plt.plot(t, dr, 'b')

    #tf = OneExpFmax()
    tf = TwoExpNoT()
    params = tf.fit_timecourse(t, dr)
    plt.plot(t, tf.fit_func(t, params), 'r')

    all_params.append(params)

all_params = np.array(all_params)

plt.figure('k')
plt.plot(np.log10(bax_concs), all_params[:,0], marker='o')
plt.figure('k linear')
plt.plot(bax_concs, all_params[:,0], marker='o')

plt.figure('fmax')
plt.plot(np.log10(bax_concs), all_params[:,1], marker='o')
plt.figure('fmax linear')
plt.plot(bax_concs, all_params[:,1], marker='o')

