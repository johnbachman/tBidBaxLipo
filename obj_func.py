#{{{# IMPORTS
from pysb import *
from pysb import anneal_sundials as anneal
from pysb.integrate import odesolve
from pylab import *
import scipy.optimize.anneal
from util.numsort import sorted_copy as sort_numeric 
from util.paramset import ParameterSet
from m0 import Vesicles_0

#from data.nbd_data import time as nbd_time, mean126c
#from grid_data import time, data_baxmaj as data
#from grid_data import time, data_baxmaj as data
import grid_analysis as g
#}}}

#{{{# GLOBALS
best_params = 0
best_err = 1e300
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
#}}}

## Annealing Functions #######################################################
#{{{# grid_anneal_fxn()
def grid_anneal_fxn(params, useparams, tmax, model, envlist, time, data, lb, ub,
                    lipo_concs_str, tbid_concs_str, bax_concs_str,
                    norm=False, vardata=False, fileobj=None):
    ''' Take the data, model, and current model parameter values at this step
        and calculate the objective function value.
    '''
    total_error = 0
    xyobs = ""
    # Check to make sure the concentrations to fit were passed in
    if (not (lipo_concs_str and tbid_concs_str and bax_concs_str)):
        raise Exception('The concentrations of tBid and Bax to fit must be specified.')

    # Check if the current parameter values are within permissible bounds
    if greater_equal(params, lb).all() and less_equal(params, ub).all():

        # Iterate over the dose response, calculating the objective function for each initial condition
        for lipo_conc_str in lipo_concs_str:
            for tbid_conc_str in tbid_concs_str:
                for bax_conc_str in bax_concs_str:
                    model.parameters['Bax_0'].value = float(bax_conc_str)
                    model.parameters['tBid_0'].value = float(tbid_conc_str)
                    model.parameters['Vesicles_0'].value = Vesicles_0.value # FIXME
     
                    # Run the model with the current initial conditions
                    outlist = anneal.annlodesolve(model, tmax, envlist, params, useparams, ic=True)
                    (xyobs, xout, yout, yobs) = outlist

                    # Calculate the objective function value
                    norm_eVes = [(eVes/Vesicles_0.value)*100 for eVes in xyobs[1]] 

                    objout = anneal.compare_data(
                               array([time, data[lipo_conc_str][tbid_conc_str][bax_conc_str]]),
                               array([xout, norm_eVes]), [(1, 1)], vardata)
                    total_error = total_error + objout
    else:
        print "======>VALUE OUT OF BOUNDS NOTED"
        temp = where((logical_and(greater_equal(params, lb), less_equal(params, ub)) * 1) == 0)
        for i in temp:
            print "======>",i, params[i]
        total_error = 1.0e300 # the largest FP in python is 1.0e308, otherwise it is just Inf

    # save the params and temps for analysis
    if fileobj and not xyobs == '':
        anneal.writetofile(fileobj, params, xyobs, objout)

    print "Total error for dose response at this step: ", total_error

    return total_error
#}}}

#{{{# nbd_anneal_fxn()
def nbd_anneal_fxn(params, useparams, tmax, model, envlist, data, lb, ub, norm=False, vardata=False, fileobj=None):
    ''' Take the nbd data, model, and current model parameter values at this step
        and calculate the objective function value.
    '''

    total_error = 0
    xyobs = ""

    norm120c = (mean120c - min(mean120c)) / (max(mean120c) - min(mean120c))

    # Check if the current parameter values are within permissible bounds
    if greater_equal(params, lb).all() and less_equal(params, ub).all():

        model.parameters['Bax_0'].value = 100 # 100 nM Bax
        model.parameters['tBid_0'].value = 20 # 20 nM tBid
        model.parameters['Vesicles_0'].value = 200000 # 200uM lipid used for NBD expt

        # Run the model with the current initial conditions
        outlist = anneal.annlodesolve(model, tmax, envlist, params, useparams, ic=True)
        (xyobs, xout, yout, yobs) = outlist
        # Calculate the objective function value

        #norm_iBax = [iBax/Vesicles_0.value)*100 for eVes in xyobs[1]] 
        iBax = xyobs[2]
        norm_iBax = (iBax - min(iBax)) / (max(iBax) - min(iBax))
        #print norm_eVes

        objout = anneal.compare_data(array([nbd_time, norm120c]), array([xout, norm_iBax]), [(1, 1)], vardata)
        total_error = objout
    else:
        print "======>VALUE OUT OF BOUNDS NOTED"
        temp = where((logical_and(greater_equal(params, lb), less_equal(params, ub)) * 1) == 0)
        for i in temp:
            print "======>",i, params[i]
        total_error = 1.0e300 # the largest FP in python is 1.0e308, otherwise it is just Inf

    # save the params and temps for analysis
    if fileobj and not xyobs == '':
        anneal.writetofile(fileobj, params, xyobs, objout)

    print "Total error for dose response at this step: ", total_error

    return total_error
#}}}

## Fitting Functions ##########################################################
#{{{# fit_grid()
def fit_grid(model, omag=3):
    data = g.data_bgsub

    #{{{# PLOT THE PRE-FIT DOSE-RESPONSE
    ### CHOOSE THE RANGE OF THE GRID TO FIT
    #lipo_concs = ['10']
    lipo_conc_str = '10' 
    lipo_concs_str = ['10']
    #tbid_concs_str = ['2', '21.5', '40']
    tbid_concs_str = ['15']

    for tbid_conc_str in tbid_concs_str:
        g.plot_timecourses(lipo_conc_str, tbid_conc_str, g.data_bgsub,
                           fixed_axis='tBid', model=model)

    bax_concs_str = g.bax_concs_str
    #}}}

    #{{{# Do the annealing =================================================
    # Initialize data structures for annealing
    # paramarr is a list of parameter values in the order generated
    # by iterating over model.parameters
    (envlist, paramarr) = anneal.annlinit(model)

    # Get the default bounds on parameters based on their nominal values
    (lb, ub, lower, upper) = anneal.getgenparambounds(paramarr, omag=omag, N=1000)

    # Run the annealing! 
    outputfile = open('anneal_output.txt', 'w')
    annlout = scipy.optimize.anneal(grid_anneal_fxn, paramarr, args=(None, 12000, model,
                        envlist, g.time, data, lb, ub, lipo_concs_str, tbid_concs_str, bax_concs_str,
                        False, False, outputfile), lower=lower, upper=upper, full_output=1)
    outputfile.close()
    #}}}

    (xmin, retval, Jmin, T, feval, iters, accept) = annlout
    # xmin is the array of best-fit parameter values. Convert these to a dict
    # indexed by parameter name:
    param_dict = {}
    for j, param in enumerate(model.parameters):
        param_dict[param.name] = xmin[j]
    # Set the retval (best error) in the annealing for this parameter set
    param_set = ParameterSet(param_dict, retval, model)

    # Plot the post-fit dose-response along with the data
    if (True):

        print "xmin: ", xmin
        print "retval: ", retval
        print "Jmin: ", Jmin
        print "T: ", T
        print "feval: ", feval
        print "iters: ", iters
        print "accept: ", accept

        for j, param in enumerate(model.parameters):
            param.value = xmin[j]

        #lipo_conc_str = '10' # For now, just use a single lipid concentration
        #tbid_conc_str = '0'
        for tbid_conc_str in tbid_concs_str:
            g.plot_timecourses(lipo_conc_str, tbid_conc_str, g.data_bgsub,
                               fixed_axis='tBid', model=model)
        
    return param_set
#}}}
