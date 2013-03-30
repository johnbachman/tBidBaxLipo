"""
Functions for normalizing, plotting, and fitting the ANTS/DPX experimental
dye release data with a variety of simple functions.

Example usage::

    import grid_analysis as g
    g.plot_timecourses('10', '21.5', g.p)

"""

import math
import numpy as np
import matplotlib.pyplot as plt
from util.fitting import Parameter, fit, fit_initial, mse
from util.numsort import sorted_copy as sort_numeric
from util.report import Report
from pysb.integrate import odesolve
from test_pandas import dt

# Select the dataset ######################
# FIXME make this a class that takes the data as an argument
print("Using gridv2.")
from data.gridv2 import time, data_tbidmaj

colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

# Assumes the grid is uniform/symmetric
# FIXME handle all this with pandas
lipo_concs_str = sort_numeric(data_tbidmaj.keys())
tbid_concs_str = sort_numeric(data_tbidmaj[lipo_concs_str[0]].keys())
bax_concs_str = sort_numeric(data_tbidmaj[lipo_concs_str[0]][tbid_concs_str[0]].keys())

# Reconfigure data array into "Bax major" order
# (i.e., indexing is lipo/Bax/tBid, rather than lipo/tBid/Bax
data_baxmaj = {}

for lipo_conc_str in lipo_concs_str:
    data_baxmaj[lipo_conc_str] = {}
    for bax_conc_str in bax_concs_str:
        data_baxmaj[lipo_conc_str][bax_conc_str] = {}
        for tbid_conc_str in tbid_concs_str:
            data_baxmaj[lipo_conc_str][bax_conc_str][tbid_conc_str] = \
                        data_tbidmaj[lipo_conc_str][tbid_conc_str][bax_conc_str]

def apply_func_to_grid(func, grid, **kw_args):
    """Apply a function to the data grid.

    Iterates through the grid structure to reach the underlying timecourse and
    apply the function, reconstructing a new data structure in which the
    timecourses have been processed by the function.
    """
    # FIXME handle this with a hierarchical dataframe, in which you iterate over the rows
    # and then apply the function to the column-based timecourse?
    out_data = {}
    x_keys = grid.keys()
    for x in x_keys:
        out_data[x] = {}
        y_keys = grid[x].keys()
        for y in y_keys:
            out_data[x][y] = {}
            z_keys = grid[x][y].keys()
            for z in z_keys:
                out_data[x][y][z] = func(grid[x][y][z], grid, x, y, z, **kw_args)
                #print out_data[x][y][z]
    return out_data

def calc_k0(timecourse, grid, x, y, z, dt=900):
    """Calculate the initial (k0) rate.

    For the given pore formation (i.e., p(t)) timecourse, calculate the initial
    pore formation rate by dividing the increase at the first time point by the
    length of the time interval (defaults to 900 seconds = 15 minutes).
    """
    # FIXME do away with this and base it on fitted functions?
    # FIXME function should take the grid and return a new dataframe,
    # FIXME where the concentration conditions are the same but instead of
    # FIXME a timecourse it returns ki, kf, and tau in the columns
    rate = (timecourse[1] - timecourse[0]) / float(dt)
    #lipo_only_timecourse = grid[x]['0']['0']
    #lipo_only_k0 = (lipo_only_timecourse[1] - lipo_only_timecourse[0]) / float(dt)
    #rate = rate - lipo_only_k0
    #bax_only_timecourse = grid[x]['0'][z]
    #bax_only_k0 = (bax_only_timecourse[1] - bax_only_timecourse[0]) / float(dt)
    #rate = rate - bax_only_k0
    return rate

def calc_ki(timecourse, grid, x, y, z, start_time_index=10):
    """Calculate the final (ki) rate.

    For the given pore formation (i.e., p(t)) timecourse, calculate a linear
    regression to the secondary kinetics, which is defined as the pore
    formation rate running from start_time_index (defaults to 1 hour, which is
    index 4) to the end.

    Returns a list containing the slope and intercept, which can also be used
    to calculate the error of the fit.
    """
    # FIXME This should go away
    tlin = time[start_time_index:len(time)]
    plin = timecourse[start_time_index:len(time)]
    (slope, intercept) = np.polyfit(tlin, plin, 1)
    return [slope, intercept]

def calc_pores(timecourse, grid, x, y, z):
    # FIXME A transformation on the data; should return another dataframe
    p_timecourse = []
    for val in timecourse:
        p  = (100 - val)/100
        if p < 0:
            p_timecourse.append(1)
            #print("E(t) > 100%!!")
        else: 
            p_timecourse.append(-math.log(float(p)))
    return p_timecourse

def calc_pores_bgsub(timecourse, grid, x, y, z):
    # FIXME A transformation on the data; should return another dataframe
    p_timecourse = []
    for i, val in enumerate(timecourse):
        # Subtract liposome only timecourse, point-by-point
        val = val - grid[x]['0']['0'][i]
        p  = (100 - val)/100
        if p < 0:
            p_timecourse.append(1)
            #p_timecourse.append(0)
            print("E(t) > 100%!!: x, y, z, i: " + x + ", " + y + ", " + z + ", " + str(i))
        elif p > 1:
            p_timecourse.append(0)
            print("E(t) < 0%!!: x, y, z, i: " + x + ", " + y + ", " + z + ", " + str(i))
        else: 
            p = -math.log(float(p))
            if (p == 0): p = 0
            p_timecourse.append(p)
            #p_timecourse.append((float(p)))
    return p_timecourse

def calc_bgsub(timecourse, grid, x, y, z):
    # FIXME A transformation on the data; should return another dataframe
    bgsub_timecourse = []
    for i, val in enumerate(timecourse):
        val = val - grid[x]['0']['0'][i]     
        bgsub_timecourse.append(val)
    return bgsub_timecourse

# FIXME Should be explicit from the setup of the dataframe
tbid_concs = [float(conc_str) for conc_str in sort_numeric(data_tbidmaj[lipo_conc_str].keys())]
bax_concs = [float(conc_str) for conc_str in sort_numeric(data_baxmaj[lipo_conc_str].keys())]

#p = apply_func_to_grid(calc_pores, data_tbidmaj)    
print("Using background (liposome-only) subtracted pore calculations.")
p = apply_func_to_grid(calc_pores_bgsub, data_tbidmaj)    
k0 = apply_func_to_grid(calc_k0, p, dt=time[1])
ki = apply_func_to_grid(calc_ki, p)
data_bgsub = apply_func_to_grid(calc_bgsub, data_tbidmaj)

## UTILITY FUNCTIONS FOR PLOTTING/ANALYSIS ###################################
# FIXME These could all go away; the caller will get the right field from the
# FIXME returned dataframe
def get_timecourses_v_bax(lipo_conc_str, tbid_conc_str):
    bax_dict = p[lipo_conc_str][tbid_conc_str]
    return [bax_dict[key] for key in sort_numeric(bax_dict.keys())]

def get_k0_v_bax(lipo_conc_str, tbid_conc_str):
    bax_dict = k0[lipo_conc_str][tbid_conc_str]
    return [bax_dict[key] for key in sort_numeric(bax_dict.keys())]

def get_k0_v_tbid(lipo_conc_str, bax_conc_str):
    bax_dict = k0[lipo_conc_str][tbid_conc_str]
    return [k0[lipo_conc_str][key][bax_conc_str] for key in tbid_concs_str]

def get_ki_v_bax(lipo_conc_str, tbid_conc_str):
    bax_dict = ki[lipo_conc_str][tbid_conc_str]
    return [bax_dict[key][0] for key in sort_numeric(bax_dict.keys())]

# FIXME This could be a good func, the timecourses could be
# FIXME calced one by one and included in a plot
def get_rate_regression(timecourse, fittype, tbid_conc=None, bax_conc=None):
    # FIXME Get rid of tbid_conc and bax_conc?
    # Initial parameter guesses
    #ki_parm = Parameter(0.0005)
    #k0_parm = Parameter(0.0015)
    #k0_parm = Parameter( (timecourse[1]-timecourse[0])/900)
    #kt_parm = Parameter(2.8e-4) # Based on a complete guess of 2500 sec for the half-life

    # Define fitting function
    #def biphasic(t): return (ki_parm()*t) + ( (k0_parm() - ki_parm()) *
    #                                          ((1 - exp(-kt_parm()*t))/kt_parm()) )
    k = Parameter(0.0025)
    k2 = Parameter(0.00025)
    fmax = Parameter(4)
    fmax2 = Parameter(0.4)
    m = Parameter(0.01)

    def single_exp (t): return ((fmax()*(1 - np.exp(-k()*t))))
    def exp_lin(t):     return ((fmax()*(1 - np.exp(-k()*t))) + (m()*t))
    def double_exp(t):  return ((fmax()*(1 - np.exp(-k()*t)))  +
                                (fmax2()*(1 - np.exp(-k2()*t))))
    def exp_exp(t):     return ((fmax()*(1 - np.exp((1- np.exp(-k()*t))   ))))

    if (fittype == 'singleexp'):
        fit(single_exp, [k, fmax], np.array(timecourse), np.array(time))
        #fit_initial(single_exp, [k, fmax], array(timecourse), array(time))
        fitfunc = single_exp
    elif (fittype == 'explin'):
        #k = Parameter(0.0033)
        #if (not (tbid_conc == None and bax_conc == None)):
        #    fmax = Parameter(0.021+(5.2e-5*(float(bax_conc)-30)))
        #    print fmax()
        fit(exp_lin, [k, fmax, m], np.array(timecourse), np.array(time))
        #fit(exp_lin, [ m], array(timecourse), array(time))
        print (fmax(), k(), m())
        #fit_initial(exp_lin, [k, fmax, m], array(timecourse), array(time))
        fitfunc = exp_lin
    elif (fittype == 'doubleexp'):
        fit(double_exp, [k, fmax, k2, fmax2], np.array(timecourse),
            np.array(time))
        #fit_initial(double_exp, [k, fmax, k2, fmax2], array(timecourse), array(time))
        fitfunc = double_exp
    elif (fittype == 'expexp'):
        fit(exp_exp, [k, fmax], np.array(timecourse), np.array(time))
        #fit_initial(exp_exp, [k, fmax], array(timecourse), array(time))
        fitfunc = exp_exp
    else:
        raise Exception('unknown fit type')

    mse_val = mse(fitfunc, np.array(timecourse), np.array(time))

    # Perform the fit
    #fit(biphasic, [ki_parm, k0_parm, kt_parm], array(timecourse), array(time))
    #fit(biphasic, [ki_parm, kt_parm], array(timecourse), array(time))
    #fit(biphasic, [ki_parm, k0_parm], array(timecourse), array(time))
    #print("k0=" + str(k0_parm()) + ", ki=" + str(ki_parm()) + ", kt=" + str(kt_parm()) )

    # Plot original values along with fitted values
    fit_time = np.linspace(0, max(time), 200)
    fit_vals = map(fitfunc, fit_time) 
    return (fit_time, fit_vals, mse_val)

## PLOTTING FUNCTIONS ########################################################
# FIXME make display optional

class GridAnalysis(object):
    """Functions for analyzing and plotting the ANTS/DPX datasets.

    Parameters
    ----------
    data : pandas.DataFrame
        A dataset from tbidbaxlipo.data.gridv1 or gridv2, in the form
        of a hierarchical dataframe. The index levels are named 'Liposomes'
        'tBid', and 'Bax' and are indexed by the concentrations/amounts
        of each, as floats. The columns of the dataframe are the timepoints,
        in seconds, as integers.
    """
    def __init__(self, data):
        self.data = data
        """The dataset to plot or analyze."""

    def plot_timecourses(self, axis_order=('tBid', 'Bax'),
                         lipo_conc=10,
                         fittype='explin', model=None,
                         report=None, display=True):
                         #major_axis=lipo_conc_str, fixed_conc_str,
                         #fixed_axis='tBid',

        if display:
            plt.ion()

        # axis_order[0] determines the concentration that we hold fixed.
        # Check what the axis is and what value we're fixing it at
        # Check which concentration axis we're varying
        # Take the cross section of the data for that concentration
        data_for_lipo_conc = self.data[lipo_conc]

        # Then get the indices for the next level
        major_axis_name = axis_order[0]
        major_axis_index = data_for_lipo_conc.columns.names.index( \
                                                            major_axis_name)
        # Iterate over the concentrations for the major axis
        for major_conc in data_for_lipo_conc.columns.levels[major_axis_index]:
            # Get a cross-section of the data for that concentration
            data_for_major_conc = data_for_lipo_conc.xs(major_conc, axis=1,
                                                        level=major_axis_name)
            # Start building a figure
            plt.figure()
            # Plot all of the timecourses across the minor axis
            color_index = 0
            for minor_conc, timecourse in data_for_major_conc.iteritems():
                plt.plot(timecourse.keys(), timecourse, 's'+colors[color_index],
                         label='%s %s' % (axis_order[1], minor_conc))
                color_index += 1
            # Add a legend
            plt.legend(loc='lower right')
            plt.title("lipo %d, %s %d" % (lipo_conc, major_axis_name,
                                          major_conc))
            # Show and/or save the figure
            if display:
                plt.show()
            if report:
                report.add_current_figure()
        # <end iteration over major-axis plots>
        # Write the report
        if report:
            report.write_report()

        """
        col_index = 0
        total_mse = 0
        if (fixed_axis == 'tBid'):
            var_concs_str = bax_concs_str
            var_axis = 'Bax'
        elif (fixed_axis == 'Bax'):
            var_concs_str = tbid_concs_str
            var_axis = 'tBid'
        else:
            raise Exception("Unknown axis: " + fixed_axis)

        # Iterate over the concentrations for the axis that varies
        for var_conc_str in var_concs_str:
            col = colors[col_index % len(colors)]
            col_index += 1
            mse_val = 0
            # Get the timecourse to plot/fit
            if (fixed_axis == 'tBid'):
                timecourse = data[lipo_conc_str][fixed_conc_str][var_conc_str]
            else:
                timecourse = data[lipo_conc_str][var_conc_str][fixed_conc_str]
            print(var_axis + ' conc ' + var_conc_str)

            # If no model object given as an argument, fit data to the specified function
            if (model == None and not fittype == None):
                (fit_time, fit_vals, mse_val) = get_rate_regression(timecourse, fittype,
                                            tbid_conc=float(fixed_conc_str),
                                            bax_conc=float(var_conc_str))
            # Otherwise, run the model with the given initial conditions
            elif (not model == None):
                fraction_dimerized = 0.02
                if (fixed_axis == 'tBid'):
                    model.parameters['Bax_0'].value = float(var_conc_str)
          #          model.parameters['Bax2_0'].value = 2 + float(fixed_conc_str)+0.001
                    model.parameters['tBid_0'].value = float(fixed_conc_str)
                else:
                    model.parameters['Bax_0'].value = float(fixed_conc_str)
          #          model.parameters['Bax2_0'].value = float(fixed_conc_str)*0 ## FIXME
                    model.parameters['tBid_0'].value = float(tBid_conc_str)

                # Run the model
                fit_time = np.linspace(0, max(time), 100)
                x = odesolve(model, fit_time)
                #fit_vals = (x['eVes']/model.parameters['Vesicles_0'].value)*100
                fit_vals = x['pores'] / 0.038
                mse_val = 0 # FIXME
            # Plot the data along with the fit/model trace
            plt.plot(time, timecourse, 's'+col, label="_nolegend_")
            if (not fittype==None):
                plt.plot(fit_time, fit_vals, '-'+col, label=var_conc_str + " " + var_axis)

            #e = abs(randn(len(timecourse)))
            #print e
            #errorbar(time, timecourse, yerr=e, fmt='s'+col)

            total_mse += mse_val

        plt.legend(loc='upper left')
        plt.title('Fits for ' + lipo_conc_str + 'uL lipid, ' + fixed_conc_str + 'nM ' + fixed_axis)
        #+ ', Fit: ' + (fittype if model == None else 'model'))
        plt.ylabel('p(t) (avg pores per vesicle)')
        plt.xlabel('Time (sec)')

        if (fittype == None):
            pass
        elif (model == None): 
            print('Fittype: ' + fittype + ', MSE: ' + str(total_mse))
        else:
            print('Fittype: model')

        if (report):
            report.addCurrentFigure()
        if display:
            plt.ioff()
    """

def test_plot_timecourses():
    """Smoke test."""
    ga = GridAnalysis(dt)
    ga.plot_timecourses(display=False, report=Report())
    assert True

def plot_dose_response(lipo_conc_str='10', rate='k0', loglogplot=False, fittype='power',
                       model=None, axis='Bax', report=None):
    """ Given a lipid concentration, plot the initial or final rate vs. Bax
    concentration dose response, with a separate curve for each tBid
    concentration.

    Parameters
    ----------
    lipo_conc_str : string
        The lipid concentration for the dose response. Default is '10'.
    rate : string: 'k0' or 'ki'
        The rate to plot as a function of concentration. Can be either
        'k0' (the initial rate) or 'ki' (the "intermediate", or final, rate).
        Default is 'k0'.
    loglogplot : boolean
        Whether to plot the dose response on a log scale or not.
        Default is False.
    fittype : None or string: 'linear', 'power', 'hill', or 'hillexp'
        The function to fit to the dose response data. If None, no fitting
        is performed. Note that fitting of this type is only done if the
        model argument is None (see below).
    model : pysb.core.Model
        The model to simulate and compare to the dose response data.
        Default is None.
    axis : string: 'tBid' or 'Bax'
        The axis over which to calculate the dose response (i.e., the y-axis).
        Default is 'Bax'.
    report : :py:class:`tbidbaxlipo.util.Report`
        The Report object to add the dose response figures to. Default is None.
    """

    plt.ion()
    plt.figure()
    col_index = 0
    total_mse = 0

    if (axis == 'Bax'):
        concs = bax_concs
        outer_concs_str = tbid_concs_str
        outer_axis = 'tBid'
    elif (axis == 'tBid'):
        concs = tbid_concs
        outer_concs_str = bax_concs_str
        outer_axis = 'Bax'
    else:
        raise Exception("Unknown axis: " + axis)

    #for tbid_conc_str in tbid_concs_str: TODO
    for outer_conc_str in outer_concs_str:
        # Choose the color for plots of this tbid concentration
        col = colors[col_index % len(colors)]
        col_index += 1

        # Check if we're plotting the initial or the secondary rate
        if (rate == 'k0'):
            if (axis == 'tBid'):
                data_arr = get_k0_v_tbid(lipo_conc_str, outer_conc_str)
            else:
                data_arr = get_k0_v_bax(lipo_conc_str, outer_conc_str)
        elif (rate == 'ki'): 
            data_arr = get_ki_v_bax(lipo_conc_str, outer_conc_str) 
        else:
            raise Exception('The rate (\'k0\' or \'ki\') was not specified.')

        model_y_vals = []
        model_x_vals = np.linspace(0, 1.01*max(concs), 50)
        fitfunc = None

        if (not model == None):
            if (axis == 'Bax'):
                for conc in model_x_vals:
                    model.parameters['tBid_0'].value = float(outer_conc_str)
                    model.parameters['Bax_0'].value = conc
                    model_y_vals.append(get_model_k0(model)) 
        # Only do fitting if there is model object passed
        else:
            # Fitting 
            if (fittype == 'linear'):
                # define fitting function
                m = Parameter(1)
                b = Parameter(0)
                def linear(x): return ((m()*x) + b())
                fit(linear, [m, b], np.array(data_arr), np.array(concs))
                #print(bid_conc_str + "nm cBid: k=" + str(k()) + ", n=" + str(n())) TODO
                fitfunc = linear
                #(slope, intercept) = np.polyfit(log_concs, log_k0, 1)
                #k0_fit = np.polyval([slope, intercept], log_concs)
                #print "slope: ", slope, ", intercept: ", intercept
                #plot(log_concs, k0_fit, '-')
            elif (fittype == 'power'):
                # define fitting function
                k = Parameter(1)
                n = Parameter(0.4)
                def powerlaw(x): return (k()*(x**n()))
                fit(powerlaw, [k, n], np.array(data_arr), np.array(concs))
                #print(bid_conc_str + "nm cBid: k=" + str(k()) + ", n=" + str(n())) TODO
                fitfunc = powerlaw
                #(slope, intercept) = np.polyfit(log_concs, log_k0, 1)
                #k0_fit = np.polyval([slope, intercept], log_concs)
                print("exponent: %d" % n())
                #plot(log_concs, k0_fit, '-')
                if (axis == 'Bax'):
                    print(outer_conc_str + "nm cBid: exponent=" + str(n()))
            elif (fittype == 'hill' or fittype == 'hillexp'):
                # define fitting function
                km = Parameter(85)
                vmax = Parameter(0.0005)
                nh = Parameter(1)
                b = data_arr[0]
                def michaelis_menten(x): return \
                        ((vmax()*(x**nh()))/((km()**nh()) + (x**nh())) + b)
                #def line(x): return (m()*x)+b()

                # Perform the fit
                if (fittype == 'hillexp'):
                    fit(michaelis_menten, [km, vmax, nh], np.array(data_arr),
                        np.array(concs))
                else:
                    fit(michaelis_menten, [km, vmax], np.array(data_arr),
                        np.array(concs))

                if (axis == 'Bax'):
                    kcat = vmax() / float(tbid_conc_str)
                    print('%s nm cBid: Km=%f, Vmax=%f, kcat=%f, nh=%f' %
                          (outer_conc_str, km(), vmax(), kcat(), nh()))
                fitfunc = michaelis_menten
            elif (fittype == None):
                pass
            else:
                raise Exception("Fitting function must be 'hill', "
                                "'hillexp', 'power', or None")

        # Plot data
        data_marker = 's-' if fittype==None else 's'
        data_legend = outer_conc_str + " " + outer_axis if fittype==None \
                                                        else '_nolegend_'

        if (loglogplot):
            plt.loglog(concs, data_arr, data_marker + col, label=data_legend)
            rise1 = np.log(data_arr[2]) - np.log(data_arr[1])
            run1 = np.log(concs[2]) - np.log(concs[1])
            slope1 = rise1 / run1
            #print "riserise
            #print "run = " + str(run)
            rise2 = np.log(data_arr[3]) - np.log(data_arr[2])
            run2 = np.log(concs[3]) - np.log(concs[2])
            slope2 = rise2 / run2
            print(outer_conc_str + 'nm ' + outer_axis + ': Slope1,2=' + str(slope1) + ', '
                            + str(slope2)) # TODO
        else:
            plt.plot(concs, data_arr, data_marker + col, label=data_legend)

        # Plot fit
        if (not model == None):
            plt.plot(model_x_vals, model_y_vals, '-'+col, label=outer_conc_str + " " + outer_axis)

        if (fitfunc):
            # Show fit error
            mse_val = mse(fitfunc, np.array(data_arr), np.array(concs))
            #print ("mse_val = " + str(mse_val))
            total_mse += mse_val
            fit_x_vals = np.linspace(0, 1.01*max(concs), 50) # TODO: Magic numbers 0 and 300
            fit_y_vals = map(fitfunc, fit_x_vals)
            if (loglogplot):
                plt.loglog(fit_x_vals, fit_y_vals, '-'+col, label=outer_conc_str + " " + outer_axis)
            else:
                plt.plot(fit_x_vals, fit_y_vals, '-'+col, label=outer_conc_str + " " + outer_axis) 

        #if (not loglogplot):
        #    plot(concs, data_arr, '-s'+col, label=tbid_conc_str + " cBid")
        # If doing log-log, ignore 0 Bax condition
        #else:
        #    log_concs = log(concs[1:len(concs)])
        #    log_data = log(data_arr[1:len(concs)])
        #    plot(log_concs, log_data, 's-'+col, label=tbid_conc_str + " cBid")

    plt.title(rate + ' vs. ' + axis + ', ' + lipo_conc_str + 'uL Lipid' +
          ('' if fittype==None else ('; Fit: ' + fittype + ', MSE: ' + str(total_mse))))
    #legend([conc + ' cBid' for conc in tbid_concs_str], loc='upper left')
    plt.legend(loc='upper left')

    plt.xlabel(axis + ' (nM)')
    plt.ylabel('rate (d(pores/ves)/dt, in seconds)')
    print("Total MSE: " + str(total_mse))

    if (report):
        report.addCurrentFigure()

def get_model_k0(model, k0_time=900):
    t = np.linspace(0, k0_time, 100)
    x = odesolve(model, t)
    return (x['pores'][-1] / model.parameters['NUM_VESICLES'].value)/k0_time



## COMMENTS ##################################################################
"""
Comments on the data: 

1. CAVEAT: CORRECTIONS TO BE MADE:

- First off, I need to REDO THE CALCULATIONS TO FIX WHERE PERMEABILIZATION
GOES OVER 100%. This should help substantially with corrections to the
p(t) calculations.

- Second, I need to subtract the permeabilization in the liposome only case.

2. OBSERVATIONS

- The time curves can be fit with two linear piecewise velocities, v0 and vi.

- The dose-response curve for k0 as a function of Bax saturates and appears to
fairly Michaelian (MAKE THESE PLOTS FOR ALL LIPID CONCS)

- The dose-response curve for ki also saturates, but in a number of instances
it actually appears to go down for the higher concentrations of Bax.
Could this indicate Bax inhibition of tBid or Bax inhibition of Bax? or tBid
inhibition of Bax? (MAKE THESE PLOTS FOR ALL LIPID CONCS)
The log-log plot of this dose response definitely appears to show at least two
limiting slopes, 1, -1, and maybe 0 also.

- For ki, the effect of cBid varies between the low lipid and high lipid conditions.
Interestingly, at high lipid concentrations, after about 2nM or 8.5nM of cBid,
all of the Bax curves lie on top of each other, indicating that after that
amount, cBid is no longer the limiting reagent.

On the other hand, for low lipid concentration, the titration lines are more
spread out between the different cBid concentrations, indicating that cBid
may be the limiting reagent, seemingly across all concentrations of Bax.

TODO List:

* Plot dose-response curves for fit values of ki, kf, t
* Plot curves for Bax v. tBid (not just tBid by Bax)
* Look at ki vs. kf, t
* Fit single funcs to all of dose response
* Average three lipid concs into one and fit
* Fit Bax only condition, for all tBids
* Get data for tBid mt1, compare

The new dataset doesn't have the vi/vf linear property?!

The new dataset has saturating velocity vs. tbid, rather than linear?

Look at tBid mt1 for both datasets.

Subtract out liposome only condition

# Things to add:
# Plot initial rate vp as a function of tbid and as a function of bax, and tbid/bax ratio
# Plot secondary rate (> 1hr?) as a function of tbid and as a function of bax, and tbid/bax ratio
"""

"""
Sources of error:
* Error in the signal read. This should be related to the actual counts generated
by the fluorimeter. For brighter (more lipid) conditions, this should be lower
to do the stronger signal.

* Error due to imprecision in the concentrations used. Was this an automated
pipettor? Presumably. Even so, there is imprecision depending on the volumes
used, especially for smaller volumes. The CV is probably 2-5%, which could make
a big difference when dealing with small differences in tBid or Bax.

* Error in the timing--probably small, though it may be large for the first time
point.

MAGNIFICATION OF ERROR
* Whether the source of error is in the initial conditions or in the signal,
the same absolute error will be magnified for high signal values. Consider
the case of 80% permeabilization vs. 15% permeabilization, with an absolute error
of 3% (permeabilization) for each.

In the case of 80% dye release, the retained dye is 20%; with 3% error, this is
17-23% retained. Calculating the avg. pores per vesicle yields:

-log(0.17) = 1.772
-log(0.23) = 1.470

With an absolute range of 0.302.

In the case of 15% dye release, the retained dye is 85%; with 3% error, this is
82-88% retained, yielding:

-log(0.82) = 0.198
-log(0.88) = 0.128

With an absolute range of 0.070.

On linear plots, the error will appear magnified for larger permeabilization values,
and for values (e.g. rates) derived from these values.

However, the magnification of error is irrelevant when fits are being performed
directly against the data.
"""

##############################################################################
### DEPRECATED/TRASH ########################################################
##############################################################################
def plot_ki_regression(lipo_conc_str, tbid_conc_str, bax_conc_str, start_time_index=4):
    """ Plot the timecourse for the given concentrations along with a linear
        regression to the latter kinetics, which is defined as the pore formation
        rate running from start_time_index (defaults to 1 hour, which is index 4)
        to the end.

        Also prints the slope of the regression and the error. """
    tlin = time[start_time_index:len(time)]
    pfull = p[lipo_conc_str][tbid_conc_str][bax_conc_str]
    plin = pfull[start_time_index:len(time)]
    (slope, intercept) = np.polyfit(tlin, plin, 1)
    pcalc = np.polyval([slope, intercept], tlin)
    plt.plot(time, pfull, 's-')
    plt.plot(tlin, pcalc)
    err = np.sqrt(sum((pcalc-plin)**2)/(len(plin)))
    print "slope: ", slope, " error: ", err
    plt.title("p(t), with k_i regression line shown")
    plt.xlabel("Time (sec)")
    plt.ylabel("p(t)")
