NBD-126-Bax Titration, Plate Reader
===================================

Data from Justin Kale, dated 9/26/2012.

Raw data
--------

By Liposome Concentration
~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::
    :context:

    from tbidbaxlipo.plots.plot_c126_titration import *
    plot_titrations_by_liposome_conc()

By Bax Concentration
~~~~~~~~~~~~~~~~~~~~

.. plot::
    :context:

    plt.close('all')
    plot_titrations_by_bax_conc()

Background-Subtracted
---------------------

.. plot::
    :context:

    plt.close('all')
    bgsub_data = background_subtract_data(raw_data)
    for bax_conc in bgsub_data.columns.levels[0]:
        plot_titration(bgsub_data, bax_conc, 'Bax')

Single-Exponential Fits
-----------------------

By Liposome Concentration
~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::
    :context:

    plt.close('all')
    k_matrix = plot_exponential_fits(bgsub_data)

By Bax and Liposome Concentration (3D Surface)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::
    :context:

    plt.close('all')
    from tbidbaxlipo.plots.plot_c126_titration import *
    plot_insertion_rate_surface(k_matrix, bgsub_data)

