Liposome Titration and Bax Insertion
====================================

Here we examine the predictions of different models on the relationship
between Bax and liposome concentration and Bax insertion kinetics.

The c126 insertion titration experiment that Justin performed seems to indicate
that increasing liposome concentrations invariably increases insertion rates,
but that increasing Bax can actually decrease these rates. It is these
characteristics that we will examine for the following models of Bax-lipid
interactions.

Partitioning model
------------------

First, we examine the partitioning model, in which the fraction of Bax
at liposomes is determined by the amount of liposomes present, and there
is no limit to how much Bax can occupy liposomes.

.. plot::

    from tbidbaxlipo.pore_comparison.schwarz import *
    from tbidbaxlipo.models import one_cpt
    plot_liposome_titration_insertion_kinetics(one_cpt)

Liposome site model
-------------------

.. plot::

    from tbidbaxlipo.pore_comparison.schwarz import *
    from tbidbaxlipo.models import lipo_sites
    plot_liposome_titration_insertion_kinetics(lipo_sites)

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

