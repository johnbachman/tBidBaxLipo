.. _stochastic_models_fig1:

Figure 1: Multi-compartment simulation of permeabilization processes
====================================================================

Figure 1A
---------
.. figure:: ../../../images/fig1a_multi_compartment.jpg
    :width: 6in
    :align: center

**Figure 1A**. In the two-compartment approach adopted by prior modeling
studies, the pore forming protein P is partitions between two reaction
compartments, solution and membrane (:math:`P_{sol}` and :math:`P_{mem}`).  The
membrane compartment represents a continuum approximation of the many discrete
vesicles in the solution. In the multi-compartment approach, individual lipid
vesicles are enumerated as explicit compartments, and the pore forming protein
partitions among each of these.

Figure 1B
---------

.. plot::
    :context:

    plt.close('all')
    from tbidbaxlipo.plots.stoch_det_comparison.bax_schwarz import jobs, data
    from tbidbaxlipo.plots.stoch_det_comparison.plots import *
    plot_dye_release_titration(jobs, data)
    xlim([0, 1500])

.. plot::
    :context:

    plt.close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pores', -1)

