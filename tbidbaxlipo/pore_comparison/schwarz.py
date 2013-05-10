import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.models.one_cpt import Builder
from pysb.integrate import odesolve

def plot_Bax_titration():

    b = Builder()
    b.build_model_bax_schwarz()

    bax_concs = np.logspace(-1, 3, 10)
    t = np.linspace(0, 12000, 1000)
    initial_rates = []

    for bax_conc in bax_concs:
        b.model.parameters['Bax_0'].value = bax_conc
        x = odesolve(b.model, t)
        avg_pores = x['pores']/b.model.parameters['Vesicles_0'].value
        initial_rate = avg_pores[50] / t[50]
        initial_rates.append(initial_rate)
    initial_rates = np.array(initial_rates)

    # Run a regression against the points and calculate the initial_rate
    log_concs = np.log(bax_concs)
    log_initial_rates = np.log(initial_rates)
    (slope, intercept) = np.polyfit(log_concs, log_initial_rates, 1)
    fitted_initial_rates = (slope * log_concs) + intercept

    # Figure with slope and intercept info
    plt.figure()
    plt.figtext(0.2, 0.8, 'Slope: %.4f' % slope)
    plt.figtext(0.2, 0.75, 'Intercept: %.4f' % intercept)
    plt.plot(log_concs, log_initial_rates, 'bx')
    plt.plot(log_concs, fitted_initial_rates, 'r')
    plt.xlabel('Log([Bax])')
    plt.ylabel('Log(V_i)')
    plt.title("Log-Log plot of initial rate vs. Bax conc")
    plt.show()

def plot_liposome_titration():
    b = Builder()
    b.translocate_Bax()
    t = np.linspace(0, 3000, 1000)
    lipo_concs = np.logspace(-2, 2, 20)
    ss_mBax_values = []

    # Plot timecourses
    plt.ion()
    plt.figure()

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
    plt.show()

if __name__ == '__main__':
    #plot_Bax_titration()
    plot_liposome_titration()
