Calculating liposome concentrations
===================================

We can roughly calculate the concentration of liposomes from the amount of
lipids we put into the reaction. For an example of the use of this type
of calculation, see [Satsoura2012]_.

First, we figure out the total amount of lipids in our solution:

.. ipython::

    In [1]: import math

    # mass of lipid film, in grams
    In [1]: mass_of_lipids = 0.002

    # resuspension volume, in liters
    In [1]: resuspension_volume = 0.001

    # concentration of lipid, in mg/ml
    In [1]: mg_per_ml_lipid = mass_of_lipids / float(resuspension_volume)

    # MW of lipids, in daltons
    In [2]: avg_mw_of_lipids = 770

    In [3]: moles_of_lipid = mass_of_lipids / float(avg_mw_of_lipids)

    In [15]: n_of_lipids = moles_of_lipid * 6.022e23

    In [16]: n_of_lipids
    Out[16]: 7.820779220779222e+17

Next, we calculate the number of lipid molecules per liposome:

.. ipython::

    # in nm
    In [5]: diameter_of_liposomes = 100

    In [9]: surface_area_of_liposomes = \
       ...: 4 * math.pi * (diameter_of_liposomes/2)**2

    # Taken from Satsoura et al. methods; in units of nm^2
    In [17]: surface_area_of_lipid_molecule = 0.75

    # We multiply by two to account for both leaflets
    In [18]: n_of_lipids_per_liposome = \
       ....: (surface_area_of_liposomes / surface_area_of_lipid_molecule) * 2

    In [19]: n_of_lipids_per_liposome
We can now calculate the number of liposomes as the number of lipid molecules
divided by the number of lipid molecules per liposome:

.. ipython::

    In [25]: n_of_liposomes = n_of_lipids / n_of_lipids_per_liposome

    In [27]: moles_of_liposomes = n_of_liposomes / 6.022e23

We are now in a position to calculate the concentration of lipids and the
concentration of liposomes given a particular resuspension volume (typically 1
mL):

.. ipython::

    In [30]: conc_of_liposomes_before_column = moles_of_liposomes / resuspension_volume

    # in Molar
    In [31]: conc_of_liposomes_before_column

    In [32]: conc_of_lipid_before_column = moles_of_lipid / resuspension_volume

    # in Molar
    In [33]: conc_of_lipid_before_column

For assays where the liposomes must be separated from free dye by running them
over a benchtop column, we must account for the dilution factor from running
the column. For a CL-2B column, we typically recover all of the liposomes
within two 1 mL fractions, for a dilution factor of 0.5 (though Aisha argues that
0.4 may be more accurate):

.. ipython::

    In [33]: dilution_factor = 0.5

    In [34]: conc_of_lipid_after_column = conc_of_lipid_before_column * dilution_factor

    In [34]: conc_of_liposomes_after_column = conc_of_liposomes_before_column * dilution_factor

    In [34]: mg_per_ml_lipid_after_column = mg_per_ml_lipid * dilution_factor

The concentration of lipids in a reaction is typically determined by the amount
of the final liposome solution is included within the 100 uL reaction volume.
Listed below is a reference table giving the final liposome concentration, in
nanomolar, for different amounts of liposome solution included in the reaction:

.. ipython::

    In [33]: from texttable import Texttable

    In [34]: tt = Texttable()

    In [35]: tt.header(['uL', '[Liposomes] (nM)', '[Lipid] (uM)',
       ....: '[Lipid] (mg/ml)'])

    In [36]: vols = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40, 50, 60, 70, \
       ....:         80, 90, 100]

    In [36]: rows = zip(vols,
       ....:  [((vol*conc_of_liposomes_after_column)/100.)*1e9 for vol in vols],
       ....:  [((vol*conc_of_lipid_after_column)/100.)*1e6 for vol in vols],
       ....:  [((vol*mg_per_ml_lipid_after_column)/100.) for vol in vols])

    In [37]: tt.add_rows(rows, header=False)

    In [38]: print tt.draw()

