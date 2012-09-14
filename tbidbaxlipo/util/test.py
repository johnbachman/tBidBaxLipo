from pysb import *

Model()

Monomer('A', ['b'])
Monomer('B', ['b'])

Parameter('k', 1)
Parameter('ic', 100)

Rule('Bind', A(b=None) + B(b=None) <> A(b=1) % B(b=1), k, k)

Initial(A(b=None), ic)
Initial(B(b=None), ic)

