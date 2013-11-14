from tbidbaxlipo.util.plate_assay import *
from tbidbaxlipo.util import error_propagation
import itertools
import pickle
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
from tbidbaxlipo.util import fitting

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict

layout = collections.OrderedDict([
        ('Bax 100 nM + 9.3 nM liposomes', ['A01', 'A03', 'A05', 'A07',
                                            'A09', 'A11']),
        ('Bax 50 nM + 9.3 nM liposomes', ['C01', 'C03', 'C05', 'C07',
                                            'C09', 'C11']),
        ('Bax 25 nM + 9.3 nM liposomes', ['E01', 'E03', 'E05', 'E07',
                                            'E09', 'E11']),
        ('Bax 12.5 nM + 9.3 nM liposomes', ['G01', 'G03', 'G05', 'G07',
                                            'G09', 'G11']),
        ('Bax 6.25 nM + 9.3 nM liposomes', ['A02', 'A04', 'A06', 'A08',
                                            'A10', 'A12']),
        ('Bax 3.13 nM + 9.3 nM liposomes', ['C02', 'C04', 'C06', 'C08',
                                            'C10', 'C12']),
        ('Bax 1.56 nM + 9.3 nM liposomes', ['E02', 'E04', 'E06', 'E08',
                                            'E10', 'E12']),
        ('Bax 0 nM + 9.3 nM liposomes', ['G02', 'G04', 'G06', 'G08',
                                            'G10', 'G12']),
        ('Bax 100 nM', ['B01', 'B03', 'B05', 'B07', 'B09', 'B11']),
        ('Bax 50 nM', ['D01', 'D03', 'D05', 'D07', 'D09', 'D11']),
        ('Bax 25 nM', ['F01', 'F03', 'F05', 'F07', 'F09', 'F11']),
        ('Bax 12.5 nM', ['H01', 'H03', 'H05', 'H07', 'H09', 'H11']),
        ('Bax 6.25 nM', ['B02', 'B04', 'B06', 'B08', 'B10', 'B12']),
        ('Bax 3.13 nM', ['D02', 'D04', 'D06', 'D08', 'D10', 'D12']),
        ('Bax 1.56 nM', ['F02', 'F04', 'F06', 'F08', 'F10', 'F12']),
        ('Bax 0 nM', ['H02', 'H04', 'H06', 'H08', 'H10', 'H12']),
    ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                 '131113_Bax488_RhoPE_37overnight.csv'))

ion()

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_wallac(timecourse_file)
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Averages of raw timecourses across replicates
(timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)
"""Averages of raw timecourses."""

# Bax only (donor-only controls)
bax_only_names = [
        'Bax 100 nM',
        'Bax 50 nM',
        'Bax 25 nM',
        'Bax 12.5 nM',
        'Bax 6.25 nM',
        'Bax 3.13 nM',
        'Bax 1.56 nM',
        'Bax 0 nM']

bax_only_wells = []
for name in bax_only_names:
    bax_only_wells += layout[name]
bax_only_set = extract(bax_only_names, timecourse_averages)

# Bax and liposomes
bax_lipos_names = [
        'Bax 100 nM + 9.3 nM liposomes',
        'Bax 50 nM + 9.3 nM liposomes',
        'Bax 25 nM + 9.3 nM liposomes',
        'Bax 12.5 nM + 9.3 nM liposomes',
        'Bax 6.25 nM + 9.3 nM liposomes',
        'Bax 3.13 nM + 9.3 nM liposomes',
        'Bax 1.56 nM + 9.3 nM liposomes',
        'Bax 0 nM + 9.3 nM liposomes',
        ]
bax_lipos_wells = []
for name in bax_lipos_names:
    bax_lipos_wells += layout[name]
bax_lipos_set = extract(bax_lipos_names, timecourse_averages)

# Liposomes only
lipos_only_avgs = extract(['Bax 0 nM + 9.3 nM liposomes'], timecourse_averages)
lipos_only_stds = extract(['Bax 0 nM + 9.3 nM liposomes'], timecourse_stds)

# Subtract liposome background (bleedthrough) from Bax+liposomes
lipo_background = np.mean(timecourse_averages['Bax 0 nM + 9.3 nM liposomes']
                                             [VALUE])
bax_lipos_bgsub = subtract_background(bax_lipos_set, lipo_background)

# Subtract buffer background from Bax only
buffer_background = np.mean(timecourse_averages['Bax 0 nM'][VALUE])
bax_bgsub = subtract_background(bax_only_set, buffer_background)


# Fret efficiency
fret_means = np.zeros(len(bax_lipos_names))
conc_list = np.zeros(len(bax_lipos_names))
fret_sds = np.zeros(len(bax_lipos_names))

for i, cond_name in enumerate(bax_lipos_names):
    conc_list[i] = float(cond_name.split(' ')[1])
    f_da = np.mean(bax_lipos_bgsub[cond_name][VALUE])
    f_da_sd = np.std(bax_lipos_bgsub[cond_name][VALUE], ddof=1)
    f_d = np.mean(bax_bgsub[bax_only_names[i]][VALUE])
    f_d_sd = np.std(bax_bgsub[bax_only_names[i]][VALUE], ddof=1)
    fret_means[i] = fret = 1 - (f_da / f_d)
    fret_sds[i] = error_propagation.calc_ratio_sd(f_da, f_da_sd, f_d, f_d_sd)

def plot_data():
    def myfig():
        figure(figsize=(9, 6))

    myfig()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    myfig()
    plot_all(timecourse_averages, errors=timecourse_stds)
    title("Raw timecourses, averaged")


if __name__ == '__main__':
    figure()
    errorbar(conc_list, fret_means, yerr=fret_sds, color='r', linewidth=2,
            linestyle='')
    title('Eq. FRET vs. [Bax]')
    xlabel('[Bax] (nM)')
    ylabel('FRET efficiency (1 - $F_{D+A}/F_D$)')
    ylim((0, 1))
    # Fit the background (donor-only) to a line
    m = fitting.Parameter(0)
    b = fitting.Parameter(0.5)
    def linear(x):
        return m()*x + b()
    fitting.fit(linear, [m, b], fret_means, conc_list)
    plot(conc_list, linear(conc_list), color='g', linewidth=2)

    figure()
    errorbar(np.log10(conc_list), fret_means, yerr=fret_sds, color='r',
            linewidth=2, linestyle='')
    title('Eq. FRET vs. [Bax]')
    xlabel('log10([Bax]) (nM)')
    ylabel('FRET efficiency (1 - $F_{D+A}/F_D$)')
    ylim((0, 1))
    plot(np.log10(conc_list), linear(conc_list), color='g', linewidth=2)

    sys.exit()

