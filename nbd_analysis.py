from data.nbd_data import *
#from data.nbd_data_62c import nbd62c, time
from util.report import Report
from util.fitting import Parameter, fit, residuals
from pylab import *


# How does signal increase with hydrophobicity? Need an additional parameter
# to describe magnitude of hydrophobicity at each position. If linear, can
# be a simple proportionality constant to the species concentration.

# Try single exponential fit to each mean curve
rep = Report()

#{{{# do_fit()
def do_fit(report=None, fittype='double_exp'):

    for i, nbd in enumerate(nbdall):
        #if (not nbd_names[i] == '3c'): continue

        if (nbd_names[i] == '62c'):
            time = time_c62
        else:
            time = time_other


        fitfig = figure()
        k1s = []
        k2s = []
        k3s = []
        fmaxs = []
        fmax2s = []
        reslist = []

        for j, replicate in enumerate(nbd):
            k1 = Parameter(0.01)
            k2 = Parameter(0.0005)
            f0 = Parameter(replicate[0])
            fmax = Parameter(0.45)
            fmax2 = Parameter(0.6)
            m = Parameter(0.01)
            k3 = Parameter(0.1)
            fmax3 = Parameter(0.0025)

            def single_exp (t):    return (f0() + (fmax()*(1 - exp(-k1()*t))))
            def exp_lin(t):        return (f0() + (fmax()*(1 - exp(-k1()*t))) + (m()*t))
            def double_exp(t):     return (f0() + (fmax()*(1 - exp(-k1()*t)))  +
                                           (fmax2()*(1 - exp(-k2()*t))))
            def triple_exp(t):     return (f0() + (fmax()*(1 - exp(-k1()*t)))  +
                                           (fmax2()*(1 - exp(-k2()*t))) +
                                           (fmax3()*(1 - exp(-k3()*t))))
            def exp_hyperbola(t):  return (f0() + (fmax()*(1 - exp(-k1()*t))) +
                                           (fmax2()*(1 - (1/(1 + (k2()*t) )) )))
            #def linked_eq(t):    return (f0() + (fmax()*exp(-k1()*t)*(-1 + exp(-k2()*t))))
            def linked_eq(t):      return (f0() + (k1()*(1 - exp(-k2()*t)) - k2()*(1 - exp(-k1()*t))) * (fmax()/(k1() - k2())))
            def linked_eq2(t):     return (f0() + (1/(k1() - k2())) * (fmax()*k1()*(exp(-k2()*t) - exp(-k1()*t)) + fmax2()*(k1()*(1 - exp(-k2()*t)) - k2()*(1 - exp(-k1()*t)))))

            def exp_exp(t):    return (f0() + (fmax()*(1 - exp(-fmax3()*(1-exp(-k3()*t))+k1()))) + (fmax2()*(1 - exp(-k2()*t))))
            #def linked_eq(t):      return (f0() + (fmax()*exp(-k3()*t)*(-1 + exp(-k2()*t)   )))

            if (fittype == 'single_exp'):
                fit(single_exp, [k1, fmax], array(replicate), array(time))
                fitfunc = single_exp
            elif (fittype == 'exp_lin'):
                fit(exp_lin, [k1, fmax, m], array(replicate), array(time))
                fitfunc = exp_lin
            elif (fittype == 'double_exp'):
                fit(double_exp, [k1, fmax, k2, fmax2], array(replicate), array(time))
                fitfunc = double_exp
            elif (fittype == 'triple_exp'):
                fit(triple_exp, [k1, fmax, k2, fmax2, k3, fmax3], array(replicate), array(time))
                fitfunc = triple_exp
            elif (fittype == 'exp_hyperbola'):
                fit(exp_hyperbola, [k1, fmax, k2, fmax2], array(replicate), array(time))
                fitfunc = exp_hyperbola
            elif (fittype == 'linked_eq'):
                fit(linked_eq, [k1, fmax, k2], array(replicate), array(time))
                fitfunc = linked_eq
            elif (fittype == 'linked_eq2'):
                fitfunc = linked_eq2
                fit(fitfunc, [k1, fmax, k2, fmax2], array(replicate), array(time))
            elif (fittype == 'exp_exp'):
                fit(exp_exp, [fmax, fmax3, k3, fmax2, k2], array(replicate), array(time))
                fitfunc = exp_exp
            else:
                raise Exception('unknown fit type')

            k1s.append(k1())
            fmaxs.append(fmax())
            k2s.append(k2())
            fmax2s.append(fmax2())
            k3s.append(k3())

            plot(time, replicate, label='No. ' + str(j), figure=fitfig)
            model_vals = map(fitfunc, time)
            plot(time, model_vals, 'k', label='__nolabel__', figure=fitfig)

            # Plot residuals
            res = residuals(fitfunc, array(replicate), array(time))
            reslist.append(res)

        legend(loc='upper right')
        xlabel('Time (seconds)')
        ylabel('Fold-Change Increase')
        title('Bax ' + nbd_names[i] + ' Data vs. ' + fittype + ' Model')

        if (report):
            report.addCurrentFigure()

        # Display some of the fitted parameter values
        print("== " + nbd_names[i] + " =============")
        print("k1:")
        print(k1s)
        print('Mean k1: %.3e, SD: %.3e' % (mean(array(k1s)), std(array(k1s))))
        print("fmax:")
        print(fmaxs)
        print('Mean fmax: %.3e, SD: %.3e' % (mean(array(fmaxs)), std(array(fmaxs))))
        print("k2:")
        print(k2s)
        print('Mean k2: %.3e, SD: %.3e' % (mean(array(k2s)), std(array(k2s))))
        print("fmax2:")
        print(fmax2s)
        print('Mean fmax2: %.3e, SD: %.3e' % (mean(array(fmax2s)), std(array(fmax2s))))

        resfig = figure()
        for res in reslist:
            plot(time, res, figure=resfig)
        title('Residuals for ' + nbd_names[i])
        if (report):
            report.addCurrentFigure()
    # end iteration over mutants

    if (report):
        report.writeReport()
#}}}#

#{{{# plot_normalized()
def plot_normalized(report=None):
    for i, nbd in enumerate(nbdall):
        if (nbd_names[i] == '62c'):
            time = time_c62
        else:
            time = time_other

        normfig = figure()
        norm_replicates = normalize(nbd)

        for j, replicate in enumerate(norm_replicates):
            plot(time, replicate, label='No. ' + str(j), figure=normfig)

        legend(loc='lower right')
        title('Bax ' + nbd_names[i] + ' Data, Normalized')

        if (report):
            report.addCurrentFigure()
    if (report):
        report.writeReport()
#}}}

#{{{# normalize()
def normalize(replicates):
    normalized_replicates = []
    for i, replicate in enumerate(replicates):
        replicate = array(replicate)
        min_val = min(replicate)
        max_val = max(replicate)
        normalized_replicate = (replicate - min_val) / (max_val - min_val)
        normalized_replicates.append(normalized_replicate)

    return normalized_replicates
#}}}

## COMMENTS #############################################################
#{{{# Fits and parameters
"""
The result of double exponential fitting yields estimates for k1 and k2
in the following order (largest/fastest to smallest/slowest):
== K1 ==
3c:   1.381e-2,  SD 9.050e-4
122c: 1.159e-2,  SD 1.140e-3
62c:  6.356e-3,  SD 1.567e-3 (nearly worthless)
126c: 5.875e-3
120c: 1.981e-3,  SD 1.233e-4
(for 126c, for some reason the first trajectory doesn't produce a good fit,
leading to a poor estimate of k1 (4.808e-3, SD 1.514e-3). A better estimate
would be to average the k1s from the two fits that actually fit:
0.00602904318 + 0.00572031140, avg: 5.875e-3

== K2 ==
3c:   1.340e-03, SD: 1.117e-04
126c: 1.139e-03 
122c: 7.454e-04, SD: 9.218e-05
120c: 1.594e-04, SD: 1.567e-04
62c:  1.098e-04, SD: 1.553e-04 (worthless)

"""
#}}}

