NBD-Bax competition assay (3/20/14)
===================================

NBD-labeled and unlabeled Bax, activated by 50 uM Bim BH3, with a fixed liposome
concentration.

.. plot::

    from tbidbaxlipo.plots.layout_140320 import *
    plt.close('all')
    plt.figure(figsize=(12, 8))
    plot_all(bgsub_wells)
    plt.title('NBD-Bax competition assay')
    plt.ylabel('NBD-Bax signal (RFU)')
