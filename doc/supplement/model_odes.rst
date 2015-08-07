Model equations
===============

.. ipython:: python

    from tbidbaxlipo.models import one_cpt
    from pysb.bng import generate_equations
    bd = one_cpt.Builder()
    bd.translocate_Bax()
    bd.basal_Bax_activation()
    generate_equations(bd.model)
    print '\n'.join(['%s:\n    %s' % tup
                     for tup in zip(bd.model.species, bd.model.odes)])
