import numpy as np
from matplotlib import pyplot as plt
from pysb.integrate import odesolve, Solver
from tbidbaxlipo.util import color_iter
from tbidbaxlipo.models import lipo_sites, one_cpt
from tbidbaxlipo.util import fitting

def plot_bax_titration(model):
    plt.ion()
    bax_concs = np.logspace(-1, 2, 40)
    t = np.linspace(0, 3000, 1000)
    initial_rates = []

    plt.figure()
    for bax_conc in bax_concs:
        model.parameters['Bax_0'].value = bax_conc
        x = odesolve(model, t)
        avg_pores = x['pores']/model.parameters['Vesicles_0'].value
        plt.plot(t, avg_pores)
        initial_rate = avg_pores[10] / t[10]
        initial_rates.append(initial_rate)
    initial_rates = np.array(initial_rates)
    plt.title('Avg. pores with Bax from 0.1 to 100')
    plt.xlabel('Time')
    plt.ylabel('Avg. pores per liposome')
    plt.show()

    # Run a regression against the points and calculate the initial_rate
    log_concs = np.log(bax_concs)
    log_initial_rates = np.log(initial_rates)
    (slope, intercept) = np.polyfit(log_concs, log_initial_rates, 1)
    fitted_initial_rates = (slope * log_concs) + intercept

    # Figure with slope and intercept info
    plt.figure()
    plt.figtext(0.2, 0.8, 'Slope: %.4f' % slope)
    plt.figtext(0.2, 0.75, 'Intercept: %.4f' % intercept)
    plt.plot(log_concs, log_initial_rates, 'b')
    plt.plot(log_concs, fitted_initial_rates, 'r')
    plt.xlabel('Log([Bax])')
    plt.ylabel('Log(V_i)')
    plt.title("Log-Log plot of initial rate vs. Bax conc")
    plt.show()

def plot_bax_titration_two_exp(model):
    plt.ion()
    bax_concs = np.logspace(-1, 3, 40)
    t = np.linspace(0, 6000, 1000)
    initial_rates = []

    fmax_list = []
    k1_list = []
    k2_list = []

    plt.figure()
    for bax_conc in bax_concs:
        model.parameters['Bax_0'].value = bax_conc
        x = odesolve(model, t)
        avg_pores = x['pores']/model.parameters['Vesicles_0'].value
        efflux = 1 - np.exp(-avg_pores)

        plt.plot(t, efflux, 'r')

        fmax = fitting.Parameter(0.9)
        k1 = fitting.Parameter(0.01)
        k2 = fitting.Parameter(0.1)

        def two_part_exp(t):
            return (fmax() * (1 - np.exp(-k1() * (1 - np.exp(-k2()*t)) * t)))

        fitting.fit(two_part_exp, [fmax, k1, k2], efflux, t)
        plt.plot(t, two_part_exp(t), 'b')

        fmax_list.append(fmax())
        k1_list.append(k1())
        k2_list.append(k2())

    fmax_list = np.array(fmax_list)
    k1_list = np.array(k1_list)
    k2_list = np.array(k2_list)

    plt.title('Avg. pores with Bax from 0.1 to 100')
    plt.xlabel('Time')
    plt.ylabel('Avg. pores per liposome')
    plt.show()

    # Run a regression against the points and calculate the initial_rate
    #log_concs = np.log(bax_concs)
    #log_initial_rates = np.log(initial_rates)
    #(slope, intercept) = np.polyfit(log_concs, log_initial_rates, 1)
    #fitted_initial_rates = (slope * log_concs) + intercept

    # k1
    plt.figure()
    plt.plot(bax_concs, k1_list, 'ro')
    plt.xlabel('[Bax] (nM)')
    plt.ylabel('$k_1$')
    plt.title("$k_1$ vs. Bax conc")
    plt.show()

    # k2
    plt.figure()
    plt.plot(bax_concs, k2_list, 'ro')
    plt.xlabel('[Bax] (nM)')
    plt.ylabel('$k_2$')
    plt.title("$k_2$ vs. Bax conc")
    plt.show()

    # Fmax 
    plt.figure()
    plt.plot(bax_concs, fmax_list, 'ro')
    plt.xlabel('[Bax] (nM)')
    plt.ylabel('$F_{max}$')
    plt.title("$F_{max}$ vs. Bax conc")
    plt.show()

def plot_liposome_titration_insertion_kinetics(module):
    """Plot the insertion kinetics of Bax over a liposome titration."""
    lipo_concs = np.logspace(-2, 2, 40)
    t = np.linspace(0, 12000, 100)

    #b_trans = module.Builder()
    #b_trans.translocate_Bax()
    b_ins = module.Builder()
    b_ins.translocate_Bax()
    b_ins.basal_Bax_activation()

    fmax_list = []
    k_list = []

    for lipo_conc in lipo_concs:
        b_ins.model.parameters['Vesicles_0'].value = lipo_conc

        # Get the SS mBax value
        #b_trans.model.parameters['Vesicles_0'].value = lipo_conc
        #s = Solver(b_trans.model, t)
        #s.run()
        #max_mBax = s.yobs['mBax'][-1]

        # Get the iBax curve
        s = Solver(b_ins.model, t)
        s.run()
        iBax = s.yobs['iBax'] / b_ins.model.parameters['Bax_0'].value
        plt.plot(t, iBax, 'r')

        # Fit to single exponential
        fmax = fitting.Parameter(0.9)
        k = fitting.Parameter(0.01)
        def single_exp(t):
            return (fmax() * (1 - np.exp(-k()*t)))
        fitting.fit(single_exp, [fmax, k], iBax, t)
        plt.plot(t, single_exp(t), 'b')

        fmax_list.append(fmax())
        k_list.append(k())

    plt.title('Inserted Bax with liposome titration')
    plt.xlabel('Time')
    plt.ylabel('Fraction of inserted Bax')
    plt.show()

    # Make plots of k and fmax as a function of lipo concentration
    fmax_list = np.array(fmax_list)
    k_list = np.array(k_list)

    # k
    plt.figure()
    plt.plot(lipo_concs, k_list, 'ro')
    plt.xlabel('Liposomes (nM)')
    plt.ylabel('$k$')
    plt.title("$k$ vs. Liposome conc")
    plt.show()

    # Fmax 
    plt.figure()
    plt.plot(lipo_concs, fmax_list, 'ro')
    plt.xlabel('Liposomes (nM)')
    plt.ylabel('$F_{max}$')
    plt.title("$F_{max}$ vs. Liposome conc")
    plt.show()

def plot_liposome_titration(builder=one_cpt.Builder()):
    b = builder
    b.translocate_Bax()
    b.model.parameters['Bax_0'].value = 100
    t = np.linspace(0, 3000, 1000)
    lipo_concs = np.logspace(-2, 4, 40)
    ss_mBax_values = []

    # Plot timecourses
    plt.ion()
    #plt.figure()

    for lipo_conc in lipo_concs:
        b.model.parameters['Vesicles_0'].value = lipo_conc
        x = odesolve(b.model, t)
        mBax_frac = x['mBax'] / b.model.parameters['Bax_0'].value
        plt.plot(t, mBax_frac)
        ss_mBax_value = mBax_frac[-1]
        ss_mBax_values.append(ss_mBax_value)
    plt.show()
    ss_mBax_values = np.array(ss_mBax_values)

    # Plot liposome titration
    plt.figure()
    plt.plot(lipo_concs, ss_mBax_values)
    plt.xlabel('Liposome concentration (nM)')
    plt.ylabel('Pct. Bax at liposomes (Bax_0 = %d nM)' %
               b.model.parameters['Bax_0'].value)
    plt.title('Bax/Liposome binding curve')
    plt.show()

def plot_fraction_bax_bound(model, figure_id=None):
    bax_concs = np.logspace(0, 4, 40)
    t = np.linspace(0, 3000, 1000)

    #plt.figure()
    ss_mBax_values = []
    for bax_conc in bax_concs:
        model.parameters['Bax_0'].value = bax_conc
        x = odesolve(model, t)
        mBax_frac = x['mBax'] / model.parameters['Bax_0'].value
        ss_mBax_value = mBax_frac[-1]
        ss_mBax_values.append(ss_mBax_value)
        #plt.plot(t, mBax_frac)
    #plt.show()

    ss_mBax_values = np.array(ss_mBax_values)

    if figure_id is not None:
        plt.figure(figure_id)
    else:
        plt.figure()

    plt.semilogx(bax_concs, ss_mBax_values, linewidth=2)
    plt.xlabel('Bax concentration (nM)')
    plt.ylabel('Pct. Bax at liposomes')
    plt.ylim([0, 1])
    plt.title('Frac Bax bound vs. Bax conc')
    plt.show()

def plot_pores_and_efflux(model):
    t = np.linspace(0, 5000, 500)
    x = odesolve(model, t)
    plt.figure()
    avg_pores = x['pores']/model.parameters['Vesicles_0'].value
    plt.plot(t, avg_pores, 'r', label="Pores/Ves")
    plt.plot(t, 1 - np.exp(-avg_pores), 'b', label="Dye release")
    plt.legend(loc='lower right')
    plt.xlabel('Time (sec)')
    plt.show()

def plot_effect_of_pore_reverse_rate():
    ci = color_iter()

    pore_reverse_rates = [1e-2, 1e-4, 1e-6]
    t = np.linspace(0, 5000, 500)
    plt.figure()

    for pore_reverse_rate in pore_reverse_rates:
        params_dict = {'Bax_0': 50., 'Vesicles_0': 50.,
                        'pore_reverse_rate_k': pore_reverse_rate}
        b = one_cpt.Builder(params_dict=params_dict)
        b.build_model_bax_schwarz_reversible()
        x = odesolve(b.model, t)
        avg_pores = x['pores']/b.model.parameters['Vesicles_0'].value
        col = ci.next()
        plt.plot(t, avg_pores, color=col, linestyle='-',
                 label="Pores, reverse rate %g" % pore_reverse_rate)
        plt.plot(t, 1 - np.exp(-avg_pores), color=col, linestyle='--',
                 label="Dye release, reverse rate %g" % pore_reverse_rate)

    plt.legend(loc='upper left')
    plt.xlabel('Time (sec)')
    plt.ylabel('Dye release/Avg. pores')
    plt.title('Dye release/pore formation with varying reverse rates')
    plt.show()

if __name__ == '__main__':
    #plot_bax_titration()
    #plot_liposome_titration()
    params_dict = {'Bax_0':100., 'Vesicles_0': 10.,
                   #'Bax_transloc_kf': 0.01,
                   #'Bax_transloc_kf': 0.1,
                   'Bax_dimerization_kf': 1e-2,
                   'Bax_dimerization_kr': 0,
                   #'Bax_tetramerization_kf': 1e-2,
                   #'Bax_tetramerization_kr': 0,
                   #'pore_formation_rate_k': 1e5,
                   #'pore_reverse_rate_k': 0
                  }

    #b = Builder()
    b = one_cpt.Builder(params_dict=params_dict)
    b.build_model_bax_heat()
    #b.build_model_bax_schwarz_tetramer_reversible()
    #plot_pores_and_efflux(b.model)
    #plot_bax_titration(b.model)
    plot_bax_titration_two_exp(b.model)
