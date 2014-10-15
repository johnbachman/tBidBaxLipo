from pysb import *

Model()

Monomer('Bax', ['conf', 'pore'], {'conf':['aq', 'mem'],
                                  'pore': ['y', 'n']})
Monomer('Vesicles', ['pore'], {'pore':['y', 'n']})

Parameter('Bax_0', 113.)
Parameter('Vesicles_0', 5.)

Parameter('Bax_to_mem_kf', 1e-3)
Parameter('Bax_to_mem_kr', 1e-1)
Parameter('Bax_forms_pores_kf', 1e-3 / Vesicles_0.value)
#Parameter('Bax_aggregates_at_pores_kf', 1e-3)

Initial(Bax(conf='aq', pore='n'), Bax_0)
Initial(Vesicles(pore='n'), Vesicles_0)

Rule('Bax_to_mem_fwd',
     Bax(conf='aq', pore='n') + Vesicles() >>
     Bax(conf='mem', pore='n') + Vesicles(),
     Bax_to_mem_kf)

Rule('Bax_to_mem_rev',
     Bax(conf='mem', pore='n') >> Bax(conf='aq', pore='n'),
     Bax_to_mem_kr)

Rule('Bax_forms_pores',
     Bax(conf='mem', pore='n') + Vesicles() >>
     Bax(conf='mem', pore='y') + Vesicles(),
     Bax_forms_pores_kf)

Rule('Vesicles_permeabilize',
     Bax(conf='mem', pore='n') + Vesicles(pore='n') >>
     Bax(conf='mem', pore='n') + Vesicles(pore='y'),
     Bax_forms_pores_kf)

"""
Rule('Bax_aggregates_at_pores',
     Bax(conf='aq', pore='n') + Vesicles(pore='y') >>
     Vesicles(pore='y'),
     Bax_aggregates_at_pores_kf)
"""

Observable('pVes', Vesicles(pore='y'))
Observable('cBax', Bax(conf='aq', pore='n'))

Expression('DR', pVes / Vesicles_0)
Expression('cBax_', cBax / Bax_0)

if __name__ == '__main__':

    from pysb.integrate import Solver
    import numpy as np
    from matplotlib import pyplot as plt
    from tbidbaxlipo.plots.stoch_det_comparison.bax_schwarz_irreversible \
            import jobs, data
    from tbidbaxlipo.plots.stoch_det_comparison.plots import *

    plt.ion()

    plot_dye_release_titration(jobs, data)

    tmax = 4000
    t = np.linspace(0, tmax, 500)
    s1 = Solver(model, t)
    s1.run()

    plt.plot(t, s1.yexpr['DR'], linewidth=2, color='g', label='pore_agg')
