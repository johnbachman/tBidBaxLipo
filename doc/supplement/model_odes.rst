Model equations
===============

.. ipython:: python

    from tbidbaxlipo.models import one_cpt
    from pysb.bng import generate_equations
    bd = one_cpt.Builder()
    bd.translocate_Bax()
    bd.basal_Bax_activation()
    generate_equations(bd.model)
    print '\n'.join(['s%d: %s:\n    %s' % (i, tup[0], tup[1])
                     for i, tup in
                        enumerate(zip(bd.model.species, bd.model.odes))])
