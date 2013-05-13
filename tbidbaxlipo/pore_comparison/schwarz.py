import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.models.one_cpt import Builder
from pysb.integrate import odesolve
from tbidbaxlipo.util import color_iter

def plot_bax_titration(model):

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

def plot_liposome_titration():
    b = Builder()
    b.translocate_Bax()
    b.model.parameters['Bax_0'].value = 10000
    t = np.linspace(0, 3000, 1000)
    lipo_concs = np.logspace(-2, 2, 40)
    ss_mBax_values = []

    # Plot timecourses
    plt.ion()
    #plt.figure()

    for lipo_conc in lipo_concs:
        b.model.parameters['Vesicles_0'].value = lipo_conc
        x = odesolve(b.model, t)
        mBax_frac = x['mBax'] / b.model.parameters['Bax_0'].value
        #plt.plot(t, mBax_frac)
        ss_mBax_value = mBax_frac[-1]
        ss_mBax_values.append(ss_mBax_value)
    #plt.show()
    ss_mBax_values = np.array(ss_mBax_values)

    # Plot liposome titration
    plt.figure()
    plt.plot(lipo_concs, ss_mBax_values)
    plt.xlabel('Liposome concentration (nM)')
    plt.ylabel('Pct. Bax at liposomes (Bax_0 = %d nM)' %
               b.model.parameters['Bax_0'].value)
    plt.title('Bax/Liposome binding curve')
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
        b = Builder(params_dict=params_dict)
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
    plt.show()

if __name__ == '__main__':
    #plot_bax_titration()
    #plot_liposome_titration()
    params_dict = {'Bax_0':100., 'Vesicles_0': 5.,
                   #'Bax_transloc_kf': 0.01,
                   #'Bax_transloc_kf': 0.1,
                   'Bax_dimerization_kf': 1e-2,
                   'Bax_dimerization_kr': 0,
                   'Bax_tetramerization_kf': 1e-2,
                   'Bax_tetramerization_kr': 0,
                   'pore_formation_rate_k': 1e5,
                   'pore_reverse_rate_k': 0
                  }

    b = Builder(params_dict=params_dict)
    b.build_model_bax_schwarz_tetramer_reversible()
    plot_pores_and_efflux(b.model)
    plot_bax_titration(b.model)
