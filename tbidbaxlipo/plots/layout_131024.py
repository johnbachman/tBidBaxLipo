from tbidbaxlipo.util.plate_assay import *
from tbidbaxlipo.util import error_propagation
import collections
import sys
import os
import tbidbaxlipo.data
import numpy as np
from tbidbaxlipo.util import fitting

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict

class Kinetics(object):
    def __init__(self, layout, timecourse_filename, num_pts_to_truncate,
                 signal_names, bg_name):
        self.layout = layout

        data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
        timecourse_file = os.path.abspath(os.path.join(data_path,
                                                       timecourse_filename))

        self.timecourse_wells = read_wallac(timecourse_file)
        self.reset_wells = reset_well_times(self.timecourse_wells)
        self.trunc_wells = truncate_timecourses(self.reset_wells,
                                                num_pts_to_truncate)

        (self.timecourse_averages, self.timecourse_stds) = \
                            averages(self.trunc_wells, self.layout)

        signal_set = extract(signal_names, self.timecourse_averages)

        # Subtract the background
        background = self.timecourse_averages[bg_name][VALUE]
        signal_wells = []
        for name in signal_names:
            signal_wells += layout[name]
        signal_wells_dict = extract(signal_wells, self.trunc_wells)
        self.bgsub_wells = subtract_background(signal_wells_dict, background)
        self.bgsub_averages = subtract_background(signal_set, background)

    def plot_raw_data(self):
        figure()
        plot_all(self.reset_wells)
        title("Bax + Liposome timecourses, raw, reset")

        figure()
        plot_all(self.trunc_wells)
        title("Bax + Liposome timecourses, first 30 points truncated")

        figure()
        plot_all(self.timecourse_averages, errors=self.timecourse_stds)

    @staticmethod
    def plot_titration(data_dict):
        conc_ss_means = []
        conc_ss_sds = []
        conc_list = []
        for conc_str in data_dict.keys():
            conc_list.append(float(conc_str.split(' ')[1]))
            conc_ss_means.append(np.mean(data_dict[conc_str][VALUE]))
            conc_ss_sds.append(np.std(data_dict[conc_str][VALUE], ddof=1))
        conc_ss_means = np.array(conc_ss_means)
        conc_ss_sds = np.array(conc_ss_sds)
        conc_list = np.array(conc_list)

        figure()
        errorbar(conc_list, conc_ss_means, yerr=conc_ss_sds, color='r')
        xlabel('[Bax] (nM)')
        ylabel('RFU')
        title('Steady-state fluorescence vs. [Bax]')

        figure()
        plot(np.log10(conc_list), np.log10(conc_ss_means), marker='o', color='r')
        xlabel('log10([Bax])')
        ylabel('RFU')
        title('Steady-state fluorescence vs. [Bax]')

        return (conc_list, conc_ss_means, conc_ss_sds)

# ==================
layout_bax_lipos = collections.OrderedDict([
        ('Bax 1089 nM + 1.55 nM liposomes', ['A01', 'B01', 'C01']),
        ('Bax 545 nM + 1.55 nM liposomes', ['A02', 'B02', 'C02']),
        ('Bax 272 nM + 1.55 nM liposomes', ['A03', 'B03', 'C03']),
        ('Bax 136 nM + 1.55 nM liposomes', ['A04', 'B04', 'C04']),
        ('Bax 68 nM + 1.55 nM liposomes', ['A05', 'B05', 'C05']),
        ('Bax 34 nM + 1.55 nM liposomes', ['A06', 'B06', 'C06']),
        ('Bax 17 nM + 1.55 nM liposomes', ['A07', 'B07', 'C07']),
        ('Bax 8.5 nM + 1.55 nM liposomes', ['A08', 'B08', 'C08']),
        ('Bax 4.3 nM + 1.55 nM liposomes', ['A09', 'B09', 'C09']),
        ('Bax 2.1 nM + 1.55 nM liposomes', ['A10', 'B10', 'C10']),
        ('Bax 1.06 nM + 1.55 nM liposomes', ['A11', 'B11', 'C11']),
        ('Bax 0 nM + 1.55 nM liposomes', ['A12', 'B12', 'C12']),
    ])
bax_lipos_signal_names = [
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
        'Bax 1.06 nM + 1.55 nM liposomes'
        ]
bax_lipos_background_name = 'Bax 0 nM + 1.55 nM liposomes'

layout_bax_only = collections.OrderedDict([
        ('Bax 1089 nM', ['D01', 'E01', 'F01']),
        ('Bax 545 nM', ['D02', 'E02', 'F02']),
        ('Bax 272 nM', ['D03', 'E03', 'F03']),
        ('Bax 136 nM', ['D04', 'E04', 'F04']),
        ('Bax 68 nM', ['D05', 'E05', 'F05']),
        ('Bax 34 nM', ['D06', 'E06', 'F06']),
        ('Bax 17 nM', ['D07', 'E07', 'F07']),
        ('Bax 8.5 nM', ['D08', 'E08', 'F08']),
        ('Bax 4.3 nM', ['D09', 'E09', 'F09']),
        ('Bax 2.1 nM', ['D10', 'E10', 'F10']),
        ('Bax 1.06 nM', ['D11', 'E11', 'F11']),
        ('Bax 0 nM', ['D12', 'E12', 'F12']),
    ])
bax_only_signal_names = [
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
        'Bax 1.06 nM'
        ]
bax_only_background_name = 'Bax 0 nM'

if __name__ == '__main__':
    ion()

    bax_only_preread = Kinetics(layout_bax_only, '131024_Bax488_preread.csv', 0,
            bax_only_signal_names, bax_only_background_name)
    bax_lipos_preread = Kinetics(layout_bax_lipos, '131024_Bax488_preread.csv', 0,
            bax_lipos_signal_names, bax_lipos_background_name)
    #preread.plot_raw_data()

    bax_lipos = Kinetics(layout_bax_lipos, '131024_Bax488_RhoPE_kinetics.csv', 30,
                         bax_lipos_signal_names, bax_lipos_background_name)
    #bax_lipos.plot_raw_data()
    (conc_list, bax_lipos_means, bax_lipos_sds) = \
                bax_lipos.plot_titration(bax_lipos.bgsub_averages)

    bax_only = Kinetics(layout_bax_only, '131024_Bax488_nolipos_kinetics.csv',
                        30, bax_only_signal_names, bax_only_background_name)
    #bax_only.plot_raw_data()
    (conc_list, bax_only_means, bax_only_sds) = \
                bax_lipos.plot_titration(bax_only.bgsub_averages)

    figure()
    plot(np.log10(conc_list), np.log10(bax_only_means), marker='o',
            color='b', label='Alexa 488 Bax')
    plot(np.log10(conc_list), np.log10(bax_lipos_means), marker='o',
            color='r', label='Alexa 488 Bax + Rho-PE liposomes')
    xlabel('log10([Bax])')
    ylabel('log10(RFU)')
    title('Steady-state fluorescence vs. [Bax]')
    legend(loc='lower right')

    # Fret by comparison with controls
    figure()
    fret = 1 - (bax_lipos_means / bax_only_means)
    num_pts = len(bax_lipos_means)
    ratio_sds = np.zeros(num_pts)
    for i in range(num_pts):
        ratio_sds[i] = error_propagation.calc_ratio_sd(
                                bax_lipos_means[i], bax_lipos_sds[i],
                                bax_only_means[i], bax_only_sds[i])
    # Log scale
    figure()
    errorbar(np.log10(conc_list), fret, yerr=ratio_sds, color='r')
    #plot(np.log10(conc_list), hill(conc_list), color='g')
    #fit_concs = 10 ** np.linspace(-1, 3.5, 100)
    #plot(np.log10(fit_concs), quad(fit_concs), color='g')
    xlabel('log10([Bax]) (nM)')
    ylabel('$F_{D+A} / F_D$')
    title('Bax/liposome FRET, 30C, 20-60sec avg')

    # -----------------------------------
    # FRET by comparison with expected F0 from fold-change
    fold_changes = []
    for well_name in bax_only.trunc_wells.keys():
        if well_name in ['D12', 'E12', 'F12']:
            continue
        preread_mean = np.mean(bax_only_preread.bgsub_wells[well_name][VALUE])
        postread_mean = np.mean(bax_only.bgsub_wells[well_name][VALUE])
        fold_change = postread_mean / preread_mean
        fold_changes.append(fold_change)
    fold_changes = np.array(fold_changes)
    figure()
    hist(fold_changes)
    mean_fold_change = np.mean(fold_changes[0:-1])
    # Calc FRET
    # Get pre_read levels, by well
    fret_fc = collections.OrderedDict()
    for well_name in bax_lipos.bgsub_wells.keys():
        preread_mean = np.mean(bax_lipos_preread.bgsub_wells[well_name][VALUE])
        exp_postread = preread_mean * mean_fold_change
        actual_bgsub_postread = np.mean(bax_lipos.bgsub_wells[well_name][VALUE])
        fret_fc[well_name] = []
        fret_fc[well_name].append(np.array([0]))
        fret_fc[well_name].append(np.array([actual_bgsub_postread/exp_postread]))
    # Calc means and standard deviations
    fret_fc_means = []
    fret_fc_sds = []
    for conc_str in bax_lipos_signal_names:
        wells = layout_bax_lipos[conc_str]
        fret_vals = [fret_fc[well_name][VALUE][0] for well_name in wells]
        fret_fc_means.append(np.mean(fret_vals))
        fret_fc_sds.append(np.std(fret_vals, ddof=1))
    fret_fc_means = np.array(fret_fc_means)
    fret_fc_sds = np.array(fret_fc_sds)
    # Fit with quadratic curve
    fmax = fitting.Parameter(1.)
    kd = fitting.Parameter(20)
    def quad(atot):
        btot = 1.55
        return 1 - fmax() * (
               (1 / (2 * atot)) *
               (atot + btot + kd() -
                np.sqrt((atot + btot + kd())**2 - (4 * atot * btot)))
               )
    fitting.fit(quad, [kd, fmax], fret_fc_means, conc_list)
    # Plot
    figure()
    errorbar(conc_list, fret_fc_means, yerr=fret_fc_sds, color='r')
    #plot(conc_list, hill(conc_list), color='g')
    plot(conc_list, quad(conc_list), color='g')
    # Log scale
    figure()
    errorbar(np.log10(conc_list), fret_fc_means, yerr=fret_fc_sds, color='r')
    #plot(np.log10(conc_list), hill(conc_list), color='g')
    fit_concs = 10 ** np.linspace(-1, 3.5, 100)
    plot(np.log10(fit_concs), quad(fit_concs), color='g')
    xlabel('log10([Bax]) (nM)')
    ylabel('$F_{D+A} / F_D$')
    title('Bax/liposome FRET, 30C, 20-60sec avg')
    sys.exit()



