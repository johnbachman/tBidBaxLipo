ANTS background in 3881 plates, bottom-read (5/22/14)
=====================================================

In the DPX requenching experiment on :ref:`5/14/14 <140514_Bax_DPX_43C>`, there
was a massively high background signal that decayed during the course of the
experiment,severely affecting the results. I inquired with Justin Kale of the
Andrews lab to see if they had observed similar things in their ANTS assays
when bottom-reading clear-bottom plates. As it turns out, they use precisely
the same plates (Corning 3881) but do not see this phenomenon. However, they
seal their plates with a clear sticker rather than a black seal, as I have been
doing. To identify the source of the background, I set up a 3881 plate with 6 wells each of water or buffer, and a black seal or no seal (uncovered). As it turns out, the presence of the seal is what gives rise to the background (see plots). In the future, I will use the clear sealing stickers that the Andrews lab uses which presumably should solve this problem.

.. plot::

    from tbidbaxlipo.plots.layout_140522 import plot_data
    plot_data()
