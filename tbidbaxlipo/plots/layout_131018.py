from tbidbaxlipo.util.plate_assay import *
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

# ~50 uL of 13.5 uM Bax + 558 uL buffer + 12 uL BSA (100x)
# yields stock conc of 1089 nM

layout = collections.OrderedDict([
        ('Bax 1089 nM + 1.55 nM liposomes', ['D01', 'E01', 'F01']),
        ('Bax 545 nM + 1.55 nM liposomes', ['D02', 'E02', 'F02']),
        ('Bax 272 nM + 1.55 nM liposomes', ['D03', 'E03', 'F03']),
        ('Bax 136 nM + 1.55 nM liposomes', ['D04', 'E04', 'F04']),
        ('Bax 68 nM + 1.55 nM liposomes', ['D05', 'E05', 'F05']),
        ('Bax 34 nM + 1.55 nM liposomes', ['D06', 'E06', 'F06']),
        ('Bax 17 nM + 1.55 nM liposomes', ['D07', 'E07', 'F07']),
        ('Bax 8.5 nM + 1.55 nM liposomes', ['D08', 'E08', 'F08']),
        ('Bax 4.3 nM + 1.55 nM liposomes', ['D09', 'E09', 'F09']),
        ('Bax 2.1 nM + 1.55 nM liposomes', ['D10', 'E10', 'F10']),
        ('Bax 1.06 nM + 1.55 nM liposomes', ['D11', 'E11', 'F11']),
        ('Bax 0 nM + 1.55 nM liposomes', ['D12', 'E12', 'F12',
                                          'G01', 'G02', 'G03', 'G04']),
        ('Bax 1089 nM', ['A01', 'B01', 'C01']),
        ('Bax 545 nM', ['A02', 'B02', 'C02']),
        ('Bax 272 nM', ['A03', 'B03', 'C03']),
        ('Bax 136 nM', ['A04', 'B04', 'C04']),
        ('Bax 68 nM', ['A05', 'B05', 'C05']),
        ('Bax 34 nM', ['A06', 'B06', 'C06']),
        ('Bax 17 nM', ['A07', 'B07', 'C07']),
        ('Bax 8.5 nM', ['A08', 'B08', 'C08']),
        ('Bax 4.3 nM', ['A09', 'B09', 'C09']),
        ('Bax 2.1 nM', ['A10', 'B10', 'C10']),
        ('Bax 1.06 nM', ['A11', 'B11', 'C11']),
        ('Bax 0 nM', ['A12', 'B12', 'C12']),
    ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                 '131018_bax.csv'))

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
        'Bax 1089 nM',
        'Bax 545 nM',
        'Bax 272 nM',
        'Bax 136 nM',
        'Bax 68 nM',
        'Bax 34 nM',
        'Bax 17 nM',
        'Bax 8.5 nM',
        'Bax 4.3 nM',
        'Bax 2.1 nM',
        'Bax 1.06 nM',
        'Bax 0 nM']
bax_only_wells = []
for name in bax_only_names:
    bax_only_wells += layout[name]
bax_only_set = extract(bax_only_names, timecourse_averages)

# Bax and liposomes
bax_lipos_names = [
        'Bax 1089 nM + 1.55 nM liposomes',
        'Bax 545 nM + 1.55 nM liposomes',
        'Bax 272 nM + 1.55 nM liposomes',
        'Bax 136 nM + 1.55 nM liposomes',
        'Bax 68 nM + 1.55 nM liposomes',
        'Bax 34 nM + 1.55 nM liposomes',
        'Bax 17 nM + 1.55 nM liposomes',
        'Bax 8.5 nM + 1.55 nM liposomes',
        'Bax 4.3 nM + 1.55 nM liposomes',
        'Bax 2.1 nM + 1.55 nM liposomes',
        'Bax 1.06 nM + 1.55 nM liposomes',
        'Bax 0 nM + 1.55 nM liposomes'
        ]
bax_lipos_wells = []
for name in bax_lipos_names:
    bax_lipos_wells += layout[name]
bax_lipos_set = extract(bax_lipos_names, timecourse_averages)

# Liposomes only
lipos_only_avgs = extract(['Bax 0 nM + 1.55 nM liposomes'], timecourse_averages)
lipos_only_stds = extract(['Bax 0 nM + 1.55 nM liposomes'], timecourse_stds)

# Subtract liposome background (bleedthrough) from Bax+liposomes
background = timecourse_averages['Bax 0 nM + 1.55 nM liposomes'][VALUE]
bax_lipos_bgsub = subtract_background(bax_lipos_set, background)

# Subtract buffer background from Bax only
buffer_background = timecourse_averages['Bax 0 nM'][VALUE]
bax_bgsub = subtract_background(bax_only_set, buffer_background)

# Fret efficiency
fret_dict = collections.OrderedDict()

for i, cond_name in enumerate(bax_lipos_names):
    time = bax_lipos_bgsub[cond_name][TIME]
    f_da = bax_lipos_bgsub[cond_name][VALUE]
    #bg = bax_only_set[bax_only_names[i]][VALUE]
    #bg = value[0]
    f_d = bax_bgsub[bax_only_names[i]][VALUE]
    fret = 1 - (f_da / f_d)
    fret_dict[cond_name] = []
    fret_dict[cond_name].append(time)
    fret_dict[cond_name].append(fret)

def plot_data():
    def myfig():
        figure(figsize=(9, 6))

    myfig()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    myfig()
    plot_all(timecourse_averages, errors=timecourse_stds)
    title("Raw timecourses, averaged")

    figure()
    lipos_only_avgs = extract(['Bax 0 nM + 1.55 nM liposomes'],
                    timecourse_averages)
    lipos_only_stds = extract(['Bax 0 nM + 1.55 nM liposomes'], timecourse_stds)
    plot_all(lipos_only_avgs, lipos_only_stds)

    myfig()
    plot_all(bax_lipos_bgsub, errors=timecourse_stds)

    myfig()
    plot_all(fret_dict, legend_position=0.7)
    ylabel('FRET efficiency')
    title('Alexa 488 L47C Bax / Rhodamine-PE liposome FRET')
"""

    myfig()
    plot_all(bgsub_timecourses, errors=timecourse_stds)
    title("Averaged and bgsub")

    myfig()
    plot_all(norm_averages, legend_position=0.7)
    title("Normalized timecourses, averaged")

    myfig()
    plot_all(pores, legend_position=0.7)
    title("Avg. pores per liposome")

    myfig()
    plot_all(treatments_cbid_bax, legend_position=0.7)
    title('Treatment comparison, cBid + Bax')

    myfig()
    plot_all(dmso, legend_position=0.7)
    title('Dye release, DMSO')

    myfig()
    plot_all(mdivi1, legend_position=0.7)
    title('Dye release, mdivi-1')

    myfig()
    plot_all(m309S, legend_position=0.7)
    title('Dye release, 309S')

    myfig()
    plot_all(m310S, legend_position=0.7)
    title('Dye release, 310S')

    myfig()
    plot_all(m365S, legend_position=0.7)
    title('Dye release, 365S')
"""

if __name__ == '__main__':
    plot_data()
    sys.exit()
    # Plot each dilution series for the Bax-only controls
    bax_only_replicates = extract(bax_only_names, layout)
    bax_only_rep_matrix = np.zeros((3, 12))
    bax_concs = []
    for i, (k, v) in enumerate(bax_only_replicates.items()):
        bax_concs.append(float(k.split(' ')[1]))
        bax_only_rep_matrix[0, i] = timecourse_wells[v[0]][VALUE][0]
        bax_only_rep_matrix[1, i] = timecourse_wells[v[1]][VALUE][0]
        bax_only_rep_matrix[2, i] = timecourse_wells[v[2]][VALUE][0]
    bax_concs = np.array(bax_concs)
    figure()
    for i in range(8):
        plot(np.log10(bax_concs), np.log10(bax_only_rep_matrix[0,:]),
             color='r', marker='o')
        plot(np.log10(bax_concs), np.log10(bax_only_rep_matrix[1,:]),
             color='g', marker='o')
        plot(np.log10(bax_concs), np.log10(bax_only_rep_matrix[2,:]),
             color='b', marker='o')

    """
    plot(np.log10(bax_concs), [np.log10(timecourse_wells[row[1]][VALUE][0])
                     for row in bax_only_replicates.values()],
                     marker='o', color='g', label='Row B')
    plot(np.log10(bax_concs), [np.log10(timecourse_wells[row[2]][VALUE][0])
                     for row in bax_only_replicates.values()],
                     marker='o', color='b', label='Row C')
    """
    legend(loc='lower right')

    # Average first 7 concentrations
    
    sys.exit()

    # Iterate over all wells in the Bax only set to get photobleaching
    # decay constant
    k_bleach_vals = []
    counter = 0
    for well in bax_only_wells:
        if counter % 3 == 0:
            figure()

        time = timecourse_wells[well][TIME]
        values = timecourse_wells[well][VALUE]
        values = values / float(values[0])
        k_bleach = fitting.Parameter(0.001)
        fmin_bleach = fitting.Parameter(0.97)
        def exp_decay(t):
            return 1 - fmin_bleach()*(1 - np.exp(-k_bleach()*t))
        fitting.fit(exp_decay, [k_bleach, fmin_bleach], values, time,
                    maxfev=1000000)
        plot(time, values, marker='o', linestyle='')
        plot(time, exp_decay(time))
        k_bleach_vals.append(k_bleach())
        counter += 1
    k_bleach_vals = np.array(k_bleach_vals)

    figure()
    hist(k_bleach_vals)

    # Fit the background (donor-only) to a line
    m = fitting.Parameter(0)
    b = fitting.Parameter(26000)
    def linear(x):
        return m()*x + b()
    t = timecourse_averages['cBid 50 nM + 0 nM liposomes'][TIME]
    x = timecourse_averages['cBid 50 nM + 0 nM liposomes'][VALUE]
    fitting.fit(linear, [m, b], x, t)
    figure()
    plot(t, x)
    plot(t, linear(t))

    # Plot each fret curve, F_DA/F_D
    figure()
    fret = collections.OrderedDict()
    f_d = bgsub_timecourses['cBid 50 nM + 0 nM liposomes'][VALUE]
    for i, conc_str in enumerate(bgsub_timecourses.keys()):
        time = bgsub_timecourses[conc_str][TIME]
        f_d_plus_a = bgsub_timecourses[conc_str][VALUE]
        f = f_d_plus_a / linear(t)
        fret[conc_str] = []
        fret[conc_str].append(time)
        fret[conc_str].append(f)
        plot(time, f, marker='o', linestyle='', color='r')

        f0 = fitting.Parameter(1)
        k = fitting.Parameter(0.001)
        fmin = fitting.Parameter(0.6)
        def exp_decay(t):
            return f0()*np.exp(-k()*t) + fmin()
        fitting.fit(exp_decay, [f0, k, fmin], f, time)
        plot(time, exp_decay(time), marker='', color='r')
    title("FRET")

    # Calculate fret at each concentration
    lipo_concs = []
    fret_0 = []
    for i, conc_str in enumerate(fret.keys()):
        #if i == 0:
        #    continue
        lipo_concs.append(float(conc_str.split(' ')[4]))
        f_0 = np.mean(fret[conc_str][VALUE][0:1])
        fret_0.append(f_0)
    lipo_concs = np.array(lipo_concs)
    fret_0 = np.array(fret_0)

    # Plot curve as down from 1
    figure()
    plot(lipo_concs, fret_0, marker='o', color='r')
    ylabel('$F_{D+A} / F_D$')
    xlabel('[Liposomes] (nM)')

    # Plot curve as up from 0
    figure()
    plot(lipo_concs, 1 - fret_0, marker='o', color='r')
    ylabel('$F_{D+A} / F_D$')
    xlabel('[Liposomes] (nM)')
    # Fit with binding curve
    kd = fitting.Parameter(3.)
    fmax = fitting.Parameter(0.25)
    def binding_eqn(x):
        return (fmax() * x) / (kd() + x)
    fitting.fit(binding_eqn, [kd, fmax], 1 - fret_0, lipo_concs)
    plot(lipo_concs, binding_eqn(lipo_concs), color='g')
    # Fit with exp curve
    kd = fitting.Parameter(1/3.)
    fmax = fitting.Parameter(0.25)
    def exp_eqn(x):
        return fmax() * (1 - np.exp(-kd()*x))
    fitting.fit(exp_eqn, [kd, fmax], 1 - fret_0, lipo_concs)
    plot(lipo_concs, binding_eqn(lipo_concs), color='m')

    # Fit with quadratic curve
    lipid_concs = lipo_concs * 83775.
    fmax = fitting.Parameter(0.25)
    kd = fitting.Parameter(3)
    def quad(atot):
        btot = 50.
        return fmax() * (
               (1 / (2 * btot)) *
               (atot + btot + kd() -
                np.sqrt((atot + btot + kd())**2 - (4 * atot * btot)))
               )
    fitting.fit(quad, [kd, fmax], 1 - fret_0, lipid_concs)
    figure()
    plot(lipid_concs, 1 - fret_0, color='r', marker='o')
    plot(lipid_concs, quad(lipid_concs), color='b')


