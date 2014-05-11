Bid/DiD liposome FRET, by Justin (4/29/14)
==========================================

1% DiD added to liposomes.

.. plot::
    :context:

    from pylab import *
    close('all')
    concs = np.array([1000, 500, 250, 125, 62.5, 31.25, 15.625,
             7.8125, 3.90625, 1.953125, 0.9765625, 0])
    fret = np.array([12.34482958, 13.32278601, 16.69284056, 18.47150343,
            19.75115244, 20.01660812, 23.22359198, 23.56273132, 27.99911577,
            29.90544125, 31.22612062, 34.06519739])
    fret_sds = np.array([1.163555487, 0.924293645, 1.928998199, 1.242210848,
                         0.762702035, 1.450931862, 1.250408076, 0.421700739,
                         0.440462041, 0.312535014, 0.700180056, 1.008105769])

    figure()
    errorbar(concs, fret, yerr=fret_sds, marker='o', color='r')
    title('Bid/DiD liposome FRET')
    xlabel('[Unlabeled Bid] (nM)')
    ylabel('FRET efficiency (%)')
    xlim([-20, 1020])
    ylim([0, 40])

.. plot::
    :context:

    close('all')
    figure()
    errorbar(log10(concs), fret, yerr=fret_sds, marker='o', color='r')
    xlabel('log10([Unlabeled Bid]) (nM)')
    ylim([0, 35])
    title('Bid/DiD liposome FRET')
