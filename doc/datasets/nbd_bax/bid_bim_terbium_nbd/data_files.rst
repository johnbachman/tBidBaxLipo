Tables of fit results
=====================

Exported data from NBD/Tb analysis.

Two conformations
-----------------

Download the fit results as
:download:`a .csv file <../../../_downloads/two_conf_fits.csv>`.

.. ipython:: python

    from tbidbaxlipo.plots.bid_bim_nbd_release import *
    import csv
    fit_results_bid = plot_2conf_fits('Bid')
    fit_results_bim = plot_2conf_fits('Bim')
    all_fit_results = fit_results_bid + fit_results_bim
    header_row = [['Activator', 'NBD Site', 'Replicate', 'Measurement',
                   'Parameter', 'Best fit value', 'Best fit std. err.']]
    csv_rows = header_row + \
                    [result.as_string_list() for result in all_fit_results]
    with open('_downloads/two_conf_fits.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(csv_rows)

Three conformations
-------------------

Download the fit results as
:download:`a .csv file <../../../_downloads/three_conf_fits.csv>`.

.. ipython:: python

    from tbidbaxlipo.plots.bid_bim_nbd_release import *
    import csv
    fit_results_bid = plot_3conf_fits('Bid')
    fit_results_bim = plot_3conf_fits('Bim')
    all_fit_results = fit_results_bid + fit_results_bim
    header_row = [['Activator', 'NBD Site', 'Replicate', 'Measurement',
                   'Parameter', 'Best fit value', 'Best fit std. err.']]
    csv_rows = header_row + \
                    [result.as_string_list() for result in all_fit_results]
    with open('_downloads/three_conf_fits.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(csv_rows)

Initial rates
-------------

Download the fit results as
:download:`a .csv file <../../../_downloads/initial_rates.csv>`.

.. ipython:: python

    from tbidbaxlipo.plots.bid_bim_nbd_release import *
    import csv
    fit_results_bid = plot_initial_rates('Bid')
    fit_results_bim = plot_initial_rates('Bim')
    all_fit_results = fit_results_bid + fit_results_bim
    header_row = [['Activator', 'NBD Site', 'Replicate', 'Measurement',
                   'Parameter', 'Best fit value', 'Best fit std. err.']]
    csv_rows = header_row + \
                    [result.as_string_list() for result in all_fit_results]
    with open('_downloads/initial_rates.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(csv_rows)







