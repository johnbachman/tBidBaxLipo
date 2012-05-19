__author__ = 'johnbachman'

from pysb import *
from pylab import *
from pysb.integrate import odesolve

t = linspace(0, 200, 100)


Model()

# Monomers
Monomer('Peptide', ['loc'], {'loc': ['w', 'lf', 'le', 'p']})
Monomer('Vesicles', ['state'], {'state': ['f', 'e']})

# Parameters
#k1_val = 1e-2       # Value assumed for pore formation rate (units of s^-1)
k1_val = 1e-1       # Value assumed for pore formation rate (units of s^-1)
v0 = 0.6/2.0
#L = 50e-6           # Total amount of lipid (in uM)
L = 10e-6          # Total amount of lipid (in uM)
beta_val = 2.6e-3   # beta = k1/k2 (unitless?)
keflx_val = 100     # Rate of dye efflux in M^-1 s^-1
P_0 = 1e-6 # initial concentration of peptide, in uM

Parameter('kon', 5.2e5) # in M^-1 s^-1
#Parameter('koff', 0.125) # in s^-1
Parameter('koff', 0) # in s^-1
Parameter('k1', k1_val)
Parameter('keflx', keflx_val/(v0*L))
Parameter('k2', k1_val / beta_val)

# Initial Conditions
Initial(Vesicles(state='f'), Parameter('Vesicles_0', L))
Initial(Peptide(loc='w'), Parameter('Peptide_0', P_0))

# Rules
Rule('Pep_Binds_fVes', Peptide(loc='w') + Vesicles(state='f') >> Peptide(loc='lf') + Vesicles(state='f'), kon)
Rule('Pep_Binds_eVes', Peptide(loc='w') + Vesicles(state='e') >> Peptide(loc='le') + Vesicles(state='e'), kon)

Rule('Pep_Unbinds_fVes', Peptide(loc='lf') >> Peptide(loc='w'), koff)
Rule('Pep_Unbinds_eVes', Peptide(loc='le') >> Peptide(loc='w'), koff)

Rule('Pep_Forms_Pores', Peptide(loc='lf') >> Peptide(loc='p'), k1)
Rule('PepPore_Releases_Dye', Peptide(loc='p') >> Peptide(loc='le'), k2)

Rule('Dye_Release', Vesicles(state='f') + Peptide(loc='p') >> Vesicles(state='e') + Peptide(loc='p'), keflx)

# Observables
Observable('eVes', Vesicles(state='e'))
Observable('pPep', Peptide(loc='p'))
Observable('lfPep', Peptide(loc='lf'))
Observable('lePep', Peptide(loc='le'))
Observable('wPep', Peptide(loc='w'))

