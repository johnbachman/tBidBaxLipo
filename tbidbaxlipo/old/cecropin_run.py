__author__ = 'johnbachman'

from pysb.integrate import odesolve
from pylab import *
from cecropin import model

def run_model():
    t = linspace(0,200,100)
    Pep_0 = 1e-6
    L_0 = 10e-6

    figure()
    x = odesolve(model,t)

    plot(t, x['eVes']/L_0 , label='eVes')
    plot(t, x['wPep']/Pep_0, label='wPep')
    plot(t, x['lfPep']/Pep_0, label='lfPep')
    plot(t, x['lePep']/Pep_0, label='lePep')
    xlabel('Time (s)')
    ylabel('Fraction of peptide species')
    #ylabel('Percentage Dye Release')
    legend()

run_model()
