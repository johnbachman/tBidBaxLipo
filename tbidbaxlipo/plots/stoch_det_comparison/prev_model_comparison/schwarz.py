from pysb import *
import math

Model()

Monomer('P', ['loc'], {'loc':['c', 'm']})
Monomer('L')
Monomer('Pore')

Parameter('P_0', 1.)
Parameter('L_0', 1.)
Parameter('P_binds_L_kf', 1.)
Parameter('P_binds_L_kr', 1.)
Parameter('Pore_formation_k', 1.)

Initial(L(), L_0)
Initial(P(loc='c'), P_0)

Rule('P_binds_L', P(loc='c') + L() >> P(loc='m') + L(), P_binds_L_kf)

Rule('P_unbinds_L', P(loc='m') >> P(loc='c'), P_binds_L_kr)

Rule('Pore_formation', P(loc='m') >> Pore(), Pore_formation_k)

Observable('Pores', Pore())

Expression('F', 1 - math.e ** (-Pores))

if __name__ == '__main__':
    from pysb.integrate import Solver
    import numpy as np
    from matplotlib import pyplot as plt

    tmax = 10
    t = np.linspace(0, tmax, 500)
    s = Solver(model, t)
    s.run()

    plt.ion()
    plt.figure()
    plt.plot(t, s.yexpr['F'])


