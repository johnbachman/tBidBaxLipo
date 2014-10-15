from pysb import *
import math

Model()

Monomer('Bax', ['conf', 'pore'], {'conf':['aq', 'mem'],
                                  'pore': ['y', 'n']})
Monomer('Vesicles', [])
Monomer('Pore', [])

Parameter('Bax_0', 100.)
Parameter('Vesicles_0', 5.)

Parameter('Bax_to_mem_kf', 5e-2)
Parameter('Bax_to_mem_kr', 1e-1)
Parameter('Bax_forms_pores_kf', 1e-3)
#Parameter('Bax_aggregates_at_pores_kf', 1e-3)

Initial(Bax(conf='aq', pore='n'), Bax_0)
Initial(Vesicles(), Vesicles_0)

Rule('Bax_to_mem_fwd',
     Bax(conf='aq', pore='n') >>
     Bax(conf='mem', pore='n'),
     Bax_to_mem_kf)

Rule('Bax_to_mem_rev',
     Bax(conf='mem', pore='n') >> Bax(conf='aq', pore='n'),
     Bax_to_mem_kr)

Rule('Bax_forms_pores',
     Bax(conf='mem', pore='n') >>
     Bax(conf='mem', pore='n') + Pore(),
     Bax_forms_pores_kf)

"""
Rule('Bax_aggregates_at_pores',
     Bax(conf='aq', pore='n') + Vesicles(pore='y') >>
     Vesicles(pore='y'),
     Bax_aggregates_at_pores_kf)
"""

Observable('pores_', Pore())

Expression('DR', 1 - math.e ** -(pores_ / Vesicles_0))

