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
from tbidbaxlipo.util.fitting import Parameter, fit, mse
from tbidbaxlipo.util.numsort import sorted_copy as sort_numeric
from tbidbaxlipo.util.report import Report
from tbidbaxlipo.data import gridv1, gridv2
from nose.tools import raises
import pandas as pd

colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'r', 'g', 'b', 'c', 'm']

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
        self.raw_data = data
        """pandas.DataFrame: The original dye release data."""
        self.bgsub_raw_data = calc_bgsub(self.raw_data)
        """pandas.DataFrame: background-subtracted raw data."""
        self.pore_data = calc_pores(self.raw_data, warnings=False)
        """pandas.DataFrame: The average number of pores per liposome."""
        self.bgsub_pore_data = calc_pores(calc_bgsub(self.raw_data))
        """pandas.DataFrame: pores calculated from bg-subtracted data."""

        self.biphasic_params = calc_biphasic_params(self.pore_data)
        """pandas.DataFrame: vi, vf, and tau for each set of concentrations."""
        self.initial_slope = calc_initial_slope(self.pore_data)
        """pandas.DataFrame: initial slope for each set of concentrations."""
        self.final_slope = calc_final_slope(self.pore_data)
        """pandas.DataFrame: final slope for each set of concentrations."""

def calc_pores(data, warnings=True):
    """Calculate the average number of pores per liposome from the data.

    This calculation, which is based on the assumption of a Poisson
    distribution of pores per liposome, can be used to calculate the
    average number of pores per liposome from the efflux signal (which
    gives the fraction of unpermeabilized liposomes). This is described
    further in the work of Gerhard Schwartz.

    Operates on the data passed in (for example, the raw data or the
    background-subtracted raw data) and returns a new dataframe containing the
    corresponding timecourses for average pore numbers.

    Parameters
    ----------
    data : pandas.DataFrame
        A dataframe containing a timecourse for each concentration, with
        the concentrations a hierarchical index on the rows and the timepoints
        as columns.
    warnings : boolean
        If True (default), warnings are printed to standard out if the
        calculated efflux value is greater than 100%.

    Returns
    -------
    pandas.DataFrame
        Hierarchical DataFrame corresponding to the original data but with
        the values in the timecourses transformed according to the formula

        pore_val = -log((100 - raw_val)/100)
    """
    # Make the pore timecourse dataframe a copy of the raw data
    pore_df = data.copy()
    # Iterate over all (liposome, tbid, bax) combinations:
    for concs, timecourse in data.iterrows():
        pore_timecourse = []
        # For every value in the raw data, convert to the -log of the
        # efflux data and append to the pore_timecourse
        for val in timecourse:
            p  = (100 - val)/100
            if p < 0:
                pore_timecourse.append(1)
                # If p < 0, then dye release is greater than 100%, so
                # print a warning
                if warnings:
                    print("Concs: %s: E(t) > 100%%!" % str(concs))
            else:
                pore_timecourse.append(-math.log(float(p)))
        # Plug the calculated pore timecourse into the pore dataframe
        pore_df.loc[concs] = pore_timecourse
    # Return the new dataframe
    return pore_df

def calc_bgsub(data):
    """Subtract the liposome-only (background) fluorescence from the data.

    Parameters
    ----------
    data : pandas.DataFrame
        A dataframe containing a timecourse for each concentration, with
        the concentrations a hierarchical index on the rows and the timepoints
        as columns.

    Returns
    -------
    pandas.DataFrame
        A dataframe corresponding to ``data``, but with the liposome-only
        timecourse subtracted from all timecourses with that liposome
        concentration.

    Raises
    ------
    ValueError
        if ``data`` does not have a index level named 'Liposomes'.
    """
    # Initialize the background-subtracted data as a copy of the original
    bgsub_data = data.copy()
    # Make sure the data has a level "Liposomes"
    if 'Liposomes' not in data.index.names:
        raise ValueError('The data does not have a level "Liposomes".')
    # Get the liposome concentrations
    liposome_axis_index = data.index.names.index('Liposomes')
    liposome_concs = data.index.levels[liposome_axis_index]

    # Iterate over the liposome concentrations
    for liposome_conc in liposome_concs:
        # Get the liposome-only timecourse
        background = data.loc[liposome_conc, 0, 0]
        # Get the cross-section of the data for this liposome concentration
        data_for_lipo_conc = data.xs(liposome_conc, axis=0, level='Liposomes')
        # Iterate over the rows in the data and subtract background
        for concs, timecourse in data_for_lipo_conc.iterrows():
            bgsub_timecourse = timecourse - background
            concs_tuple = tuple([liposome_conc] + list(concs))
            bgsub_data.loc[concs_tuple] = bgsub_timecourse

    return bgsub_data

def calc_initial_slope(data):
    """Calculate the initial (vi) rate.

    For the given pore formation (i.e., p(t)) timecourse, calculate the
    initial pore formation rate by dividing the increase at the first time
    point by the length of the time interval.

    Parameters
    ----------
    data : pandas.DataFrame
        A dataframe containing a timecourse for each concentration, with
        the concentrations a hierarchical index on the rows and the timepoints
        as columns.

    Returns
    -------
    pandas.DataFrame
        A dataframe containing an initial slope value for each concentration
        combination. The dataframe has the single column named 'vi',
        while the concentrations are a hierarchical index on the rows.
    """
    # Initialize a dict to hold the values
    initial_slope_dict = {}
    # Iterate over all (liposome, tbid, bax) combinations:
    for concs, timecourse in data.iterrows():
        # Calculate the initial slope and put it in the dict
        timepoints = timecourse.keys().values
        initial_rate = ((timecourse.values[1] - timecourse.values[0]) /
                        float(timepoints[1] - timepoints[0]))
        initial_slope_dict[concs] = pd.Series({'vi': initial_rate})
    df = pd.DataFrame(initial_slope_dict).T
    df.index = data.index
    return df

def calc_final_slope(data, start_time_index=10):
    """Calculate the final (vf) rate.

    For the given pore formation (i.e., p(t)) timecourse, calculate a linear
    regression to the secondary kinetics, which is defined as the pore
    formation rate running from start_time_index (defaults to 10,
    which regresses over the last four timepoints).

    Parameters
    ----------
    data : pandas.DataFrame
        A dataframe containing a timecourse for each concentration, with
        the concentrations a hierarchical index on the rows and the timepoints
        as columns.

    Returns
    -------
    pandas.DataFrame
        A dataframe containing a final slope value and intercept for each
        concentration combination. The dataframe has two columns named 'vf'
        and 'intercept', while the concentrations are a hierarchical index
        on the rows.
    """
    # Initialize a dict to hold the values
    final_slope_dict = {}
    # Iterate over all (liposome, tbid, bax) combinations:
    for concs, timecourse in data.iterrows():
        # Calculate the final slope and put it in the dict
        timepoints = timecourse.keys().values
        # Run a linear regression
        tlin = timepoints[start_time_index:len(timepoints)]
        plin = timecourse.values[start_time_index:len(timecourse.values)]
        (slope, intercept) = np.polyfit(tlin, plin, 1)
        final_slope_dict[concs] = pd.Series({'vf': slope,
                                             'intercept': intercept})
    df = pd.DataFrame(final_slope_dict).T
    df.index = data.index
    return df

def calc_biphasic_params(data):
    r"""Calculate a two-phase exponential fit for each set of concentrations.

    See Schwarz et al., "Kinetics of melittin induced pore formation in the
    membranes of lipid vesicles," Biochimica et Biophysica Acta, 1110 (1992)
    97-104.

    Fits the equation

    .. math::

        v_f t + (v_i - v_f) \left(\frac{1 - e^{-\tau t}}{\tau} \right)

    to the data.

    Parameters
    ----------
    data : pandas.DataFrame
        A dataframe containing a timecourse for each concentration, with
        the concentrations a hierarchical index on the rows and the timepoints
        as columns.

    Returns
    -------
    pandas.DataFrame
        Hierarchical DataFrame corresponding to the original data but with
        the data for "vi", "vf", and "tau" in the columns.
    """
    # A dict to hold the parameters
    biphasic_dict = {}
    # Iterate over all (liposome, tbid, bax) combinations:
    for concs, timecourse in data.iterrows():
        # Get an exponential linear fit and put the parameters
        # in the dict
        biphasic_fit_result = get_timecourse_fit(timecourse,
                                                 fittype='biphasic')
        biphasic_dict[concs] = pd.Series(biphasic_fit_result.parameters)
    # Initialize a DataFrame from the dict
    df = pd.DataFrame(biphasic_dict).T
    df.index = data.index
    return df

class FitResult(object):
    """A class to store the results of fitting a model to a timecourse.

    Parameters
    ----------
    fit_time : numpy.array
        An array of timepoints for the fitted function.
    fit_vals : numpy.array
        The predicted values at each timepoint.
    mse_val : number
        Mean squared error of the fit.
    parameters : dict
        Parameters of the fit, as name/value pairs.
    """
    def __init__(self, fit_time, fit_vals, mse_val, parameters):
        self.fit_time = fit_time
        self.fit_vals = fit_vals
        self.mse_val = mse_val
        self.parameters = parameters

def get_timecourse_fit(timecourse, fittype='biphasic'):
    """Get a fit curve for a timecourse.

    Parameters
    ----------
    timecourse : pandas.Series
        The timecourse to be fit.
    fittype : string
        One of the following:
            * `singleexp`. Single exponential.
            * `biphasic`. Two-phase exponential (see Schwarz).
            * `explin`. Exponential plus linear term (default).
            * `doubleexp`. Sum of two exponentials
            * `expexp`. Exponential raised to an exponential.

    Returns
    -------
    FitResult object
        Contains the (time, val) coordinates for the fitted curve as well
        as the mean squared error and parameter values.
    """

    # Initial parameter guesses
    k = Parameter(0.0025)
    k2 = Parameter(0.00025)
    fmax = Parameter(4)
    fmax2 = Parameter(0.4)
    m = Parameter(0.01)

    #vi = Parameter( (timecourse.values[1]-timecourse.values[0])/900)
    vi = Parameter(0.0005)
    vf = Parameter(0.0005)
    # Based on a complete guess of 2500 sec for the half-life
    tau = Parameter(2.8e-4)

    # Define fitting function
    def biphasic(t):    return (vf()*t) + ( (vi() - vf()) *
                                            ((1 - np.exp(-tau()*t))/tau()) )
    def single_exp(t):  return ((fmax()*(1 - np.exp(-k()*t))))
    def exp_lin(t):     return ((fmax()*(1 - np.exp(-k()*t))) + (m()*t))
    def double_exp(t):  return ((fmax()*(1 - np.exp(-k()*t)))  +
                                (fmax2()*(1 - np.exp(-k2()*t))))
    def exp_exp(t):     return ((fmax()*(1 - np.exp((1- np.exp(-k()*t))   ))))

    parameters = None
    # Run the fit
    if (fittype == 'biphasic'):
        fit(biphasic, [vi, vf, tau], timecourse.values,
            timecourse.keys().values)
        fitfunc = biphasic
        parameters = {'vi': vi(), 'vf': vf(), 'tau': tau()}
    elif (fittype == 'singleexp'):
        fit(single_exp, [k, fmax], timecourse.values, timecourse.keys().values)
        fitfunc = single_exp
        parameters = {'k': k(), 'fmax':fmax()}
    elif (fittype == 'explin'):
        fit(exp_lin, [k, fmax, m], timecourse.values, timecourse.keys().values)
        fitfunc = exp_lin
        parameters = {'k': k(), 'fmax':fmax(), 'm':m()}
    elif (fittype == 'doubleexp'):
        fit(double_exp, [k, fmax, k2, fmax2],
            timecourse.values, timecourse.keys().values)
        fitfunc = double_exp
        parameters = {'k': k(), 'fmax':fmax(), 'k2':k2(), 'fmax':fmax2()}
    elif (fittype == 'expexp'):
        fit(exp_exp, [k, fmax], timecourse.values, timecourse.keys().values)
        fitfunc = exp_exp
        parameters = {'k': k(), 'fmax':fmax()}
    else:
        raise Exception('unknown fit type')

    # Calculate the mean squared error of the fit
    mse_val = mse(fitfunc, timecourse.values, timecourse.keys().values)

    # Return time/value pairs for fit curve, along with the error
    fit_time = np.linspace(0, max(timecourse.keys().values), 200)
    fit_vals = map(fitfunc, fit_time) 
    return FitResult(fit_time, fit_vals, mse_val, parameters)

def plot_timecourses(data,
                     axis_order=('Liposomes', 'tBid', 'Bax'),
                     fixed_conc=10,
                     fittype='biphasic', model=None,
                     report=None, display=True):
    """Plot a cross-section of the ANTS/DPX timecourses.

    In generating the plots, one axis (the "fixed" axis), is presumed
    to be held constant, (e.g., the amount of liposomes). The next
    axis is the "major" axis: each value along the major axis (e.g.,
    a concentration of tBid) is used to generate a figure. The values
    along the "minor" axis (e.g., concentration of Bax) give rise to
    curves in each figure.

    Which concentration (liposomes, tBid, Bax) should serve as the
    fixed, major, and minor axes can be specified by the arguments to the
    function.

    Parameters
    ----------
    data : pandas.DataFrame
        A dataframe with timecourses for each (Liposome, tBid, Bax) 
        concentration combination. The data could be the raw data,
        average pore per liposome data, background subtracted, etc.
        The concentrations should be a hierarchical index on the rows,
        with the time points as the columns.
    axis_order : 3-tuple of strings
        The order of axes. axis_order[0] denotes the fixed axis,
        axis_order[1] the major axis, and axis_order[2] the minor axis.
    fixed_conc : number
        The concentration index into the fixed axis (specified by
        ``axis_order[0]``) for the value to use for all plots (e.g., 10,
        for the 10 uL condition for the liposome axis).
    fittype : string
        A string (to be passed to :py:func:`get_timecourse_fit`)
        specifying the curve fit type to be plot along with the data.
    model : pysb.core.Model or None
        If a mechanistic model is to be fit to the data, this can be set
        here. If present it overrides the fittype argument.
    report : tbidbaxlipo.util.Report
        If not None, the report to save the figures to. Default is None.
    display : boolean
        Specifies whether the figures should be displayed interactively.
        Default is True.

    Raises
    ------
    ValueError
        * if the fixed_conc specified is not a valid key for the fixed
          axis (i.e., axis_order[0])
        * if one of the axis names in axis_order is not found in the
          levels of the row index.
    """

    # Make sure the fixed_conc specified is a valid key for the fixed axis
    fixed_axis_index = data.index.names.index(axis_order[0])
    data.index.levels[fixed_axis_index]
    if fixed_conc not in data.index.levels[fixed_axis_index]:
        raise ValueError('fixed_conc %s is not a valid index into the '
                         'axis level %s!' % (fixed_conc, axis_order[0]))
    # Make sure all of the axes are legitimate
    for axis_name in axis_order:
        if axis_name not in data.index.names:
            raise ValueError('axis %s not in the index of the levels of the '
                             'data.' % axis_name)

    # If the plots should be displayed, turn on interactive mode
    if display:
        plt.ion()

    # axis_order[0] determines the concentration that we hold fixed.
    # Take the cross section of the data for the fixed concentration:
    data_for_fixed_conc = data.xs(fixed_conc, axis=0, level=axis_order[0])

    # Then get the level index for the next level
    major_axis_name = axis_order[1]
    major_axis_index = data_for_fixed_conc.index.names.index( \
                                                        major_axis_name)
    # Iterate over the concentrations for the major axis
    for major_conc in data_for_fixed_conc.index.levels[major_axis_index]:

        # Get a cross-section of the data for the major concentration
        data_for_major_conc = data_for_fixed_conc.xs(major_conc, axis=0,
                                                     level=major_axis_name)

        # -- Plot all of the timecourses across the minor axis --
        plt.figure()
        color_index = 0
        for minor_conc, timecourse in data_for_major_conc.iterrows():
            # Plot the data points
            plt.plot(timecourse.keys().values, timecourse,
                     's'+colors[color_index], label='__nolegend__')
            # -- Fit data --
            # If no model object given as an argument, fit data to the
            # specified function
            if (model == None and not fittype == None):
                fit_result = get_timecourse_fit(timecourse, fittype)
                plt.plot(fit_result.fit_time, fit_result.fit_vals,
                         '-'+colors[color_index],
                         label="%s %s" % (axis_order[2], minor_conc))
            color_index += 1
        # Add a legend
        plt.legend(loc='lower right')
        plt.title("Fits for %.2f uL liposomes, %d nM %s " % \
                    (fixed_conc, major_conc, major_axis_name))
                #+ ', Fit: ' + (fittype if model == None else 'model'))
        #plt.ylabel('p(t) (avg pores per vesicle)') # FIXME
        plt.xlabel('Time (sec)')
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
        # Otherwise, run the model with the given initial conditions
        elif (not model == None):
            fraction_dimerized = 0.02
            if (fixed_axis == 'tBid'):
                model.parameters['Bax_0'].value = float(var_conc_str)
      #          model.parameters['Bax2_0'].value = 2 + float(fixed_conc_str)+0.001
                model.parameters['tBid_0'].value = float(fixed_conc_str)
            else:
                model.parameters['Bax_0'].value = float(fixed_conc_str)
      #          model.parameters['Bax2_0'].value = float(fixed_conc_str)*0
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

def plot_titration(data, param_name, axis_order=('Liposomes', 'tBid', 'Bax'),
                   display=True, report=None):
    """Plot kinetic parameters vs. concentration of tBid, Bax, or liposomes.

    Plots the value of a kinetic parameter as a function of concentration.
    Takes a DataFrame with parameter names (from a fit) as the columns,
    and the concentrations as hierarchical indices on the rows. For each
    concentration in the major axis, (e.g., liposomes), makes a figure; for
    each concentration in the minor axis (e.g. tBid), adds a curve to the
    figure; the value of the parameter is then plotted with the innermost
    axis (e.g., Bax) as the x-coordinate.

    Parameters
    ----------
    data : pandas.DataFrame
        The parameter data derived from a fit to the timecourse data.
    param_name : string
        The name of the parameter to plot vs. concentration.
    axis_order : 3-tuple of strings
        The order of axes. axis_order[0] denotes the major axis,
        axis_order[1] the minor axis, and axis_order[2] the innermost
        (x-coordinate) axis.
    report : tbidbaxlipo.util.Report
        If not None, the report to save the figures to. Default is None.
    display : boolean
        Specifies whether the figures should be displayed interactively.
        Default is True.

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

    Raises
    ------
    ValueError
        * if the param_name specified is not a valid column name for the
          given parameter data.
        * if one of the axis names in axis_order is not found in the
          levels of the row index.
    """

    # Make sure the named parameter is a valid column in the data table
    if param_name not in data.columns:
        raise ValueError('Parameter named %s is not a valid column in '
                         'the data table.' % param_name)
    # Make sure all of the axes are legitimate
    for axis_name in axis_order:
        if axis_name not in data.index.names:
            raise ValueError('axis %s not in the index of the levels of the '
                             'data.' % axis_name)
    # If the plots should be displayed, turn on interactive mode
    if display:
        plt.ion()

    # Get the major axis
    major_axis_name = axis_order[0]
    major_axis_index = data.index.names.index(major_axis_name)
    major_axis_concs = data.index.levels[major_axis_index]
    # Get the minor axis
    minor_axis_name = axis_order[1]
    minor_axis_index = data.index.names.index(minor_axis_name)
    minor_axis_concs = data.index.levels[minor_axis_index]

    total_mse = 0

    # Iterate over the concentrations for the major axis
    for major_conc in major_axis_concs:
        plt.figure()
        color_index = 0
        # Iterate over the concentrations for the minor axis
        for minor_conc in minor_axis_concs:
            # Get the double cross section for the major/minor combination
            titration_params = data.xs((major_conc, minor_conc), axis=0,
                                   level=(major_axis_name, minor_axis_name))
            titration = titration_params[param_name]

            # Plot the data points
            plt.plot(titration.keys().values, titration,
                     's-'+colors[color_index])

            # -- Fit titration --
            # If no model object given as an argument, fit data to the
            # specified function
            """
            if (model == None and not fittype == None):
                fit_result = get_timecourse_fit(timecourse, fittype)
                plt.plot(fit_result.fit_time, fit_result.fit_vals,
                         '-'+colors[color_index],
                         label="%s %s" % (axis_order[2], minor_conc))
            """
            color_index += 1
        # Add a legend to the figure
        #plt.legend(loc='lower right')
        plt.title("Titrations for %f %s " % (major_conc, major_axis_name))
                #+ ', Fit: ' + (fittype if model == None else 'model'))
        #plt.ylabel('p(t) (avg pores per vesicle)') # FIXME
        #plt.xlabel('Time (sec)')
        #plt.title(rate + ' vs. ' + axis + ', ' + lipo_conc_str + 'uL Lipid' +
        #      ('' if fittype==None else ('; Fit: ' + fittype +
        #      ', MSE: ' + str(total_mse))))
        #legend([conc + ' cBid' for conc in tbid_concs_str], loc='upper left')
        #plt.legend(loc='upper left')
        #plt.xlabel(axis + ' (nM)')
        #plt.ylabel('rate (d(pores/ves)/dt, in seconds)')
        #print("Total MSE: " + str(total_mse))

        # Show and/or save the figure
        if display:
            plt.show()
        if report:
            report.add_current_figure()
    # <end iteration over major-axis plots>
    # Write the report
    if report:
        report.write_report()

def get_titration_fit():
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
            #print(bid_conc_str + "nm cBid: k=" + str(k()) + ", n=" + str(n())) 
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
            #print(bid_conc_str + "nm cBid: k=" + str(k()) + ", n=" + str(n())) 
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
            print(outer_conc_str + 'nm ' + outer_axis +
                    ': Slope1,2=' + str(slope1) + ', ' + str(slope2)) # TODO
        else:
            plt.plot(concs, data_arr, data_marker + col, label=data_legend)

        # Plot fit
        if (not model == None):
            plt.plot(model_x_vals, model_y_vals, '-'+col,
                    label=outer_conc_str + " " + outer_axis)

        if (fitfunc):
            # Show fit error
            mse_val = mse(fitfunc, np.array(data_arr), np.array(concs))
            #print ("mse_val = " + str(mse_val))
            total_mse += mse_val
            # TODO: Magic numbers 0 and 300
            fit_x_vals = np.linspace(0, 1.01*max(concs), 50)
            fit_y_vals = map(fitfunc, fit_x_vals)
            if (loglogplot):
                plt.loglog(fit_x_vals, fit_y_vals, '-'+col,
                        label=outer_conc_str + " " + outer_axis)
            else:
                plt.plot(fit_x_vals, fit_y_vals, '-'+col,
                            label=outer_conc_str + " " + outer_axis) 

        #if (not loglogplot):
        #    plot(concs, data_arr, '-s'+col, label=tbid_conc_str + " cBid")
        # If doing log-log, ignore 0 Bax condition
        #else:
        #    log_concs = log(concs[1:len(concs)])
        #    log_data = log(data_arr[1:len(concs)])
        #    plot(log_concs, log_data, 's-'+col, label=tbid_conc_str + " cBid")

## TESTS ###########

def test_calc_pores():
    """GridAnalysis.calc_pores should run without error."""
    calc_pores(gridv1.data, warnings=False)
    assert True

def test_calc_bgsub():
    """GridAnalysis.calc_bgsub should run without error."""
    calc_bgsub(gridv1.data)
    assert True

@raises(ValueError)
def test_calc_bgsub_no_liposomes():
    """Should raise an error if the data does not have a level 'Liposomes'."""
    calc_bgsub(gridv1.data.loc[10])

def test_calc_initial_slope():
    """GridAnalysis.calc_initial_slope() should run without error."""
    calc_initial_slope(gridv1.data)
    assert True

def test_calc_final_slope():
    """GridAnalysis.calc_final_slope() should run without error."""
    calc_final_slope(gridv1.data)
    assert True

def test_calc_explin_params():
    """Should run without error."""
    calc_biphasic_params(calc_pores(gridv1.data))

def test_get_timecourse_fit():
    """GridAnalysis.get_timecourse_fit should run without error."""
    get_timecourse_fit(gridv1.data.loc[10,0,0], fittype='singleexp')
    get_timecourse_fit(gridv1.data.loc[10,0,0], fittype='explin')
    get_timecourse_fit(gridv1.data.loc[10,0,0], fittype='doubleexp')
    get_timecourse_fit(gridv1.data.loc[10,0,0], fittype='expexp')
    assert True

def test_plot_timecourses():
    """plot_timecourses should run without error."""
    ga = GridAnalysis(gridv1.data)
    plot_timecourses(ga.raw_data, display=False)
    plot_timecourses(ga.pore_data, display=False)
    assert True

@raises(ValueError)
def test_plot_timecourses_illegal_fixed_conc():
    """plot_timecourses should error if the fixed_conc is an invalid key
    for the fixed axis (i.e., for axis_order[0])."""
    ga = GridAnalysis(gridv1.data)
    plot_timecourses(ga.raw_data, display=False,
                     axis_order=('tBid', 'Liposomes', 'Bax'), fixed_conc=10)

@raises(ValueError)
def test_plot_timecourses_illegal_level_name():
    """plot_timecourses should error if one of the level names in axis_order
    is not present in the DataFrame."""
    ga = GridAnalysis(gridv1.data)
    plot_timecourses(ga.raw_data, display=False,
                     axis_order=('tBid', 'Liposomes', 'bad_axis_name'),
                     fixed_conc=10)

def test_plot_titration():
    """plot_titration should run without error."""
    ga = GridAnalysis(gridv1.data)
    plot_titration(ga.biphasic_params, 'vi')
    assert True

@raises(ValueError)
def test_plot_titration_illegal_param_name():
    """plot_titration should error if the param_name is an invalid column
    name in the fit_params DataFrame."""
    ga = GridAnalysis(gridv1.data)
    plot_titration(ga.biphasic_params, 'bad_name', display=False,
                   axis_order=('tBid', 'Liposomes', 'Bax'))

@raises(ValueError)
def test_plot_titration_illegal_level_name():
    """plot_titration should error if one of the level names in axis_order
    is not present in the DataFrame."""
    ga = GridAnalysis(gridv1.data)
    plot_titration(ga.biphasic_params, 'vi', display=False,
                   axis_order=('tBid', 'Liposomes', 'bad axis name'))


# TODO could add test for raising an exception for an unknown fit type

# -- REFACTORING LINE --




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
