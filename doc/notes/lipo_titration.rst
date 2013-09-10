Liposome titration and Bax insertion
====================================

Here we examine the predictions of different models on the relationship
between Bax and liposome concentration and Bax insertion kinetics.

The c126 insertion titration experiment that Justin performed seems to indicate
that increasing liposome concentrations invariably increases insertion rates,
but that increasing Bax can actually decrease these rates. It is these
characteristics that we will examine for the following models of Bax-lipid
interactions.

.. _bax-insertion-partitioning-model:

Partitioning model
------------------

**Conclusion**: TOO SIMPLISTIC

First, we examine the partitioning model, in which the fraction of Bax at
liposomes is determined by the amount of liposomes present, and there is no
limit to how much Bax can occupy liposomes.

.. plot::

    from tbidbaxlipo.plots.pore_plots import *
    from tbidbaxlipo.models import one_cpt
    b = one_cpt.Builder()
    b.translocate_Bax()
    b.basal_Bax_activation()
    plot_liposome_titration_insertion_kinetics(b.model)

In this model, increasing the amount of Bax has no affect on the insertion
kinetics at all:

.. plot::

    from tbidbaxlipo.plots.pore_plots import *
    from tbidbaxlipo.models import one_cpt
    b = one_cpt.Builder()
    b.translocate_Bax()
    b.basal_Bax_activation()
    plot_bax_titration_insertion_kinetics(b.model)

Liposome site model
-------------------

**Conclusion:** IMPLAUSIBLE

It appears that in the liposome binding site model, inserted Bax can plateau at
a level well bellow 100% of Bax, as liposome binding sites are occupied by
irreversibly inserted Bax, preventing further Bax from binding and inserting.

So...the results from the c126 Bax titration assay seem to look more like
the partition model than the site model, insofar as low amounts of liposomes
relative to Bax appear not to result in low steady state levels of insertion.

One question: the data here is fold-change data over background. Since the "background" level likely incorporates fluorescence from fully aqueous as well as
peripheral forms of Bax, with likely differences in fluorescence, the
background may not scale linearly with the total amount of Bax (though
this could be measured fairly in a no-tBid titration experiment).

.. plot::

    from tbidbaxlipo.plots.pore_plots import *
    from tbidbaxlipo.models import lipo_sites
    b = lipo_sites.Builder()
    b.translocate_Bax()
    b.basal_Bax_activation()
    plot_liposome_titration_insertion_kinetics(b.model)

For a titration of Bax, the liposome site model predicts that the Fmax for
insertion should decrease sharply with increased Bax, because a smaller and
smaller fraction of total Bax gets bound and activated.

.. plot::

    from tbidbaxlipo.plots.pore_plots import *
    from tbidbaxlipo.models import lipo_sites
    b = lipo_sites.Builder()
    b.translocate_Bax()
    b.basal_Bax_activation()
    plot_bax_titration_insertion_kinetics(b.model)

.. _bax-insertion-tbid-activation:

tBid Activates Bax
------------------

**Conclusion**: PRETTY GOOD

A model with tBid-activation of Bax has the right shape for the liposome
titration:


.. plot::

    from tbidbaxlipo.plots.pore_plots import *
    from tbidbaxlipo.models import one_cpt
    params_dict = {'tBid_mBax_kf':1e-1, 'tBid_iBax_kc':1e-2}
    b = one_cpt.Builder(params_dict)
    b.translocate_tBid_Bax()
    b.tBid_activates_Bax()
    plot_liposome_titration_insertion_kinetics(b.model)

Interestingly, setting the parameters for tBid activation to have a low
Km can reproduce the phenomenon of increasing Bax causing slower kinetics:

.. plot::

    from tbidbaxlipo.plots.pore_plots import *
    from tbidbaxlipo.models import one_cpt
    params_dict = {'tBid_mBax_kf':1e-1, 'tBid_iBax_kc':1e-2}
    b = one_cpt.Builder(params_dict)
    b.translocate_tBid_Bax()
    b.tBid_activates_Bax()
    plot_bax_titration_insertion_kinetics(b.model)

tBid Activates and Binds Bax
----------------------------

**Conclusion**: DUBIOUS

Adding Bax inhibition of tBid (by tBid-activated Bax binding) produces
titration curves that are faster with increasing liposome concentration, as
occurs in the c126 titration data. However, the individual kinetic curves have
the characteristic late-linear slope due to inhibition of the enzyme, and
the curves from the plate reader do not have this characteristic (it is
worth noting that the original data from the plate reader does appear to
have a bit of this behavior, however).

.. plot::

    from tbidbaxlipo.plots.pore_plots import *
    from tbidbaxlipo.models import one_cpt
    b = one_cpt.Builder()
    b.translocate_tBid_Bax()
    b.tBid_activates_Bax()
    b.iBax_binds_tBid_at_bh3()
    plot_liposome_titration_insertion_kinetics(b.model)

Similarly, the Bax titration shows that the insertion kinetics get slower with
increasing Bax, but again, the curves have a late-linear slope that does not
appear in the data. Moreover, the fitted Fmax values go down sharply, which
should be a testable prediction.

.. plot::

    from tbidbaxlipo.plots.pore_plots import *
    from tbidbaxlipo.models import one_cpt
    b = one_cpt.Builder()
    b.translocate_tBid_Bax()
    b.tBid_activates_Bax()
    b.iBax_binds_tBid_at_bh3()
    plot_bax_titration_insertion_kinetics(b.model)

