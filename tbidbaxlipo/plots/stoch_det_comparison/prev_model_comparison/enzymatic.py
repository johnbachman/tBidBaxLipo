from pysb import *

Model()

Monomer('P')
Monomer('L', ['dye'], {'dye':['f', 'e']})

Parameter('kobs', 1.)
Parameter('P_0', 1.)
Parameter('L_0', 1.)

Initial(P(), P_0)
Initial(L(dye='f'), L_0)

Rule('P_permeabilizes_L',
     P() + L(dye='f') >> P() + L(dye='e'),
     kobs)

Observable('L_empty', L(dye='e'))

Expression('F', L_empty / L_0)

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

