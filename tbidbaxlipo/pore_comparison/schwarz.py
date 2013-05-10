import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.models.one_cpt import Builder
from pysb.integrate import odesolve

bax_concs = np.logspace(-1, 3, 10)
b = Builder()
b.build_model_bax_schwarz()

def plot_Bax_titration():
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
    plt.figure(2)
    plt.figtext(0.2, 0.8, 'Slope: %.4f' % slope)
    plt.figtext(0.2, 0.75, 'Intercept: %.4f' % intercept)
    plt.plot(log_concs, log_initial_rates, 'bx')
    plt.plot(log_concs, fitted_initial_rates, 'r')
    plt.xlabel('Log([Bax])')
    plt.ylabel('Log(V_i)')
    plt.title("Log-Log plot of initial rate vs. Bax conc")
    plt.show()

if __name__ == '__main__':
    plot_Bax_titration()
