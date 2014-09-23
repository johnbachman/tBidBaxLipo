from pysb import *

Model()

Monomer('L', ['dye'], {'dye':['f', 'e']})
Monomer('P', ['loc', 'dye', 'pore'],
        {'loc':['c', 'm'], 'dye':['none', 'f','e'], 'pore':['y', 'n']})

Parameter('P_L_kf', 1.)
Parameter('P_L_kr', 1.)
Parameter('P_forms_pores_k', 1.)
Parameter('Pores_close_k', 1.)
Parameter('Dye_release_k', 1.)
Parameter('L_0', 1.)
Parameter('P_0', 1.)

Initial(P(loc='c', dye='none', pore='n'), P_0)
Initial(L(dye='f'), L_0)

Rule('P_binds_Lf',
     P(loc='c', dye='none') + L(dye='f') >> P(loc='m', dye='f') + L(dye='f'),
     P_L_kf)

Rule('P_binds_Le',
     P(loc='c', dye='none') + L(dye='e') >> P(loc='m', dye='e') + L(dye='e'),
     P_L_kf)

Rule('P_unbinds_L',
     P(loc='m', pore='n') >> P(loc='c', pore='n', dye='none'),
     P_L_kr)

Rule('P_forms_pores',
     P(loc='m', pore='n', dye='f') >> P(loc='m', pore='y', dye='f'),
     P_forms_pores_k)

Rule('Pores_close',
     P(loc='m', pore='y', dye='f') >> P(loc='m', pore='n', dye='e'),
     Pores_close_k)

Rule('Dye_release',
     L(dye='f') + P(pore='y') >> L(dye='e') + P(pore='y'),
     Dye_release_k)

Observable('L_empty', L(dye='e'))

Expression('F', L_empty / L_0)

if __name__ == '__main__':
    from pysb.integrate import Solver
    import numpy as np
    from matplotlib import pyplot as plt

    tmax = 1000
    t = np.linspace(0, tmax, 500)
    s = Solver(model, t)
    s.run()

    plt.ion()
    plt.figure()
    plt.plot(t, s.yexpr['F'])


