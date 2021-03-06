"""
A set of functions useful for analyzing and plotting data from timecourse
assays from the Wallac plate reader.

Plate-based kinetics
--------------------

- Export the data, sorted by well. To do this, take the Microsoft Excel output
  file, go to the `List` tab; highlight the data; then go to the `Data` menu
  and select `Sort by Well`. Finally export the data to CSV as a `Windows
  Comma Separated (.csv)` file.

.. note:: .csv Compatibility

    The csv reader module used in this script appears to have a problem with
    `.csv` files exported from Excel using the default `.csv` settings.
    However, it works fine with files exported using the `Windows .csv`
    setting.

- Create a plate layout dict, for example, the following::

    expt = {'No tBid, No Bax' : ['A01', 'B01'], # Leaving out C01
            'No tBid, Bax' : ['D01', 'E01', 'F01'],
            'tBid 1, No Bax' : ['A02', 'B02', 'C02'],
            'tBid 1, Bax' : ['D02', 'E02', 'F02'],
            'tBid 2, No Bax' : ['A03', 'B03', 'C03'],
            'tBid 2, Bax' : ['D03', 'E03', 'F03'],
            'tBid 3, No Bax' : ['A04', 'B04'], # Leaving out C04
            'tBid 3, Bax' : ['D04', 'E04', 'F04'],
            'Bim BH3, No Bax' : ['A05', 'B05', 'C05'],
            'Bim BH3, Bax' : ['D05', 'E05', 'F05'],
            'Cecropin A' : ['G01', 'G02', 'G03']}

Well-based kinetics
-------------------

For well-based kinetics assays, the data does not need to be sorted by well,
since all of the measurements for the well appear in a set before the
measurements for the next well. Instead, the data can simply be read with
:py:func:`read_wallac` after saving the Excel file to `.csv`. However, note
that the instrument runs a continuous timer for all kinetics readings, so
if the kinetics reading is initiated by dispensing a reagent (e.g., liposomes)
then the time coordinates will need to be explicitly reset for each well.
This can be done using the function :py:func:`reset_well_times`.

"""
# TODO Write a function that will read the initial and final .csv or
# excel files and use the values to normalize the timecourses

import csv
import datetime as dt
from matplotlib.pyplot import *
from matplotlib.font_manager import FontProperties
import numpy as np
import collections
import pandas as pd
from copy import copy

# We store the well data in a list of lists; the first list
# has the time coordinates, the second list has the values
# (see read_wallac, below)
TIME = 0
VALUE = 1

def read_envision(csv_file):
    return read_wallac(csv_file, WELL_COL=4, TIME_COL=10, VALUE_COL=13)

def read_wallac(csv_file, WELL_COL=2, TIME_COL=4, VALUE_COL=5):
    """Parses the plate reader data from a CSV file.

    NOTE: The data must be sorted by well before exporting to CSV. See the
    instructions in the module docstring for more details.

    Parameters
    ----------
    csv_file : string
        The name of the CSV file to open.

    Returns
    -------
    dict
        A dict mapping the name of the well to a list of lists. For each dict
        entry, the first element in the outer list is a list of the time
        coordinates associated with each timepoint; the second element in the
        outer list is the list of assay values for each timepoint.  For
        example, ``wells['A01']`` returns::

            [[79.8, 145.43, ..., 2639.35, 2704.98],
             [2285,   2318, ...,    2323,    2304]]
    """

    # The plan: iterate over the CSV file, collecting time coordinates
    # and values until we find that we're at a new well; at which point
    # we update the value of cur_well, make a new entry in the dict,
    # and continue.
    csv_reader = csv.reader(open(csv_file, 'r'))

    # Initialize the empty dict
    wells = collections.OrderedDict([])
    # We're not currently on any well!
    cur_well = None
    # Iterate over the rows:
    for row in csv_reader:
        well_name = row[WELL_COL]

        # If we're at the header line, skip it
        if well_name == 'Well':
           continue
        # If it's a new well, or if the well entry is blank (as it is in
        # kinetics measurements, with multiple readings for the well), then
        # update cur_well and initialize the dict with a pair of empty lists
        elif well_name != '' and cur_well != well_name:
            cur_well = well_name
            wells[well_name] = [np.array([]),np.array([])]

        # Add the time/val pair to the dict
        time = dt.datetime.strptime(row[TIME_COL], '%H:%M:%S.%f')
        secs = (time.hour * 3600) + (time.minute * 60) + time.second + \
               (time.microsecond * 0.000001)
        value = int(row[VALUE_COL])
        wells[cur_well][TIME] = np.append(wells[cur_well][TIME], secs)
        wells[cur_well][VALUE] = np.append(wells[cur_well][VALUE], value)

    return wells

def read_flexstation_flex(tsv_file):
    """Reads a Flex mode datafile exported from the FlexStation 3.

    NOTE: Data should be organized by time, so that each row contains
    a timepoint and each column contains data for one well.

    In Flex mode files, time values for the reads are contained in distinct
    columns with names like "A8T" with time values for reads in column "A8".
    To parse these we first treat these time columns as wells like the others,
    then at the end we rebuild the dict and put the corresponding time and
    value columns together.

    Parameters
    ----------
    tsv_file : string
        Path to the tab-separated file to load.

    Returns
    -------
    collections.OrderedDict
        An ordered dictionary with well names as keys and two-element lists
        as the entries. The first entry in each two-element list is a list
        of time coordinates; the second is a list of values (e.g.,
        fluorescence values).
    """
    # The plan: get the names of all the wells by iterating over the header
    # column, then populate the entries for those wells by iterating
    csv_reader = csv.reader(open(tsv_file), dialect='excel-tab')
    # Initialize the empty dict
    time_array = []
    raw_wells = collections.OrderedDict([])
    for row_counter, row in enumerate(csv_reader):
        # Check and see if we've finished reading all the data (which would
        # show up as a blank line/empty row)
        if not row:
            break
        # Skip the first two rows
        if row_counter < 2:
            continue
        # If we're in the second row, which contains the well names, create
        # entries in the dict for the wells
        elif row_counter == 2:
            # Save the header row so that we can use it for indexing later
            header_row = row
            # Now iterate over each column in the row
            for i, well_name in enumerate(row):
                # Skip the first two columns, which contain the plate and
                # temperature headings, respectively
                if i < 2:
                    continue
                # Create an empty list to which we will add the time/value lists
                else:
                    raw_wells[well_name] = []
        # For the rows containing data, fill in the data values; the time arrays
        # will be added later
        else:
            for i, value in enumerate(row):
                if i < 2:
                    continue
                # Add the time or fluorescence value to the appropriate well,
                # unless it's empty, in which case we add nothing
                else:
                    if value == '':
                        continue
                    header_name = header_row[i]
                    if value == '#Sat':
                        raw_wells[header_name].append(np.nan)
                    elif value == '#Low':
                        raw_wells[header_name].append(np.nan)
                    else:
                        raw_wells[header_name].append(float(value))
    # Now that we're done iterating over the file, iterate over the dict we've
    # created and add the time array to each entry in the dict
    wells = collections.OrderedDict([])
    for well_name, values in raw_wells.iteritems():
        # If this well doesn't contain any data, skip it
        if not values:
            continue
        # If we're in a time well, parse out the well name
        if well_name.endswith('T'):
            index = TIME
            well_name = well_name[:-1]
        else:
            index = VALUE

        # If this is a new well, create it
        if not well_name in wells.keys():
            wells[well_name] = [[],[]]
        # Add the data to the appropriate field for the well
        wells[well_name][index] = np.array(values)
    return wells

def read_flexstation_kinetics(tsv_file):
    """Reads a kinetics datafile exported from the FlexStation 3.

    NOTE: Data should be organized by time, so that each row contains
    a timepoint and each column contains data for one well.

    Parameters
    ----------
    tsv_file : string
        Path to the tab-separated file to load.

    Returns
    -------
    collections.OrderedDict
        An ordered dictionary with well names as keys and two-element lists
        as the entries. The first entry in each two-element list is a list
        of time coordinates; the second is a list of values (e.g.,
        fluorescence values).
    """
    # The plan: get the names of all the wells by iterating over the header
    # column, then populate the entries for those wells by iterating
    csv_reader = csv.reader(open(tsv_file), dialect='excel-tab')
    # Initialize the empty dict
    time_array = []
    wells = collections.OrderedDict([])
    for row_counter, row in enumerate(csv_reader):
        # Check and see if we've finished reading all the data (which would
        # show up as a blank line/empty row)
        if not row:
            break
        # Skip the first two rows
        if row_counter < 2:
            continue
        # If we're in the second row, which contains the well names, create
        # entries in the dict for the wells
        elif row_counter == 2:
            # Save the header row so that we can use it for indexing later
            header_row = row
            # Now iterate over each column in the row
            for i, well_name in enumerate(row):
                # Skip the first two columns, which contain the time and
                # temperature headings, respectively
                if i < 2:
                    continue
                # Create an empty list to which we will add the time/value lists
                else:
                    wells[well_name] = [np.array([]), np.array([])]
        # For the rows containing data, fill in the data values; the time arrays
        # will be added later
        else:
            for i, value in enumerate(row):
                # Add the time value to the time array, after formatting
                if i == 0:
                    # If there is one colon (':') separator, then there is no hour
                    # prefix
                    if value.count(':') == 1:
                        time = dt.datetime.strptime(value, '%M:%S')
                    # If there are two colon separators, there is an hour prefix
                    elif value.count(':') == 2:
                        time = dt.datetime.strptime(value, '%H:%M:%S')
                    # Otherwise, we're confused
                    else:
                        raise Exception('Could not parse time field: too many colons!')

                    secs = (time.hour * 3600) + (time.minute * 60) + \
                            time.second + (time.microsecond * 0.000001)
                    time_array = np.append(time_array, secs)
                # Skip the temperature column
                elif i == 1:
                    continue
                # Otherwise, add the value to the appropriate well, unless
                # it's empty, in which case we add nothing
                else:
                    if value == '':
                        continue
                    well_name = header_row[i]
                    if value == '#Sat':
                        wells[well_name][VALUE] = \
                                np.append(wells[well_name][VALUE], np.nan)
                    else:
                        wells[well_name][VALUE] = \
                                np.append(wells[well_name][VALUE], float(value))

    # Now that we're done iterating over the file, iterate over the dict we've
    # created and add the time array to each entry in the dict
    for well_name, entry in wells.iteritems():
        entry[TIME] = copy(time_array)
        if len(entry[VALUE]) == 0:
            del wells[well_name]
    return wells

def plot_all(wells, errors=None, do_legend=True, legend_position=0.8):
    """Plots all of the timecourses in the dict.

    Parameters
    ----------
    wells : dict of timecourse data of the type returned by read_wallac.
    """
    for wellname, wellval in wells.iteritems():
        if errors is None:
            plot(wellval[TIME], wellval[VALUE], label=wellname)
        else:
            errorbar(wellval[TIME], wellval[VALUE],
                    yerr=errors[wellname][VALUE], label=wellname)

    # Label the axes and add a legend
    xlabel('Time') # TODO: automatically determine seconds or minutes, etc.
    ylabel('Value')

    if do_legend:
        fontP = FontProperties()
        fontP.set_size('small')
        ax = gca()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * legend_position,
                         box.height])
        legend(loc='upper left', prop=fontP, ncol=1, bbox_to_anchor=(1, 1),
             fancybox=True, shadow=True)


def plot_subset(wells, well_name_list):
    """Plots the timecourses associated with a subset of wells.

    Parameters
    ----------
    wells : dict
        dict of timecourse data of the type returned by read_wallac.
    well_name_list : list of strings
        The names of the wells to plot.
    """

    for well_name in well_name_list:
        timecourse = wells[well_name]
        plot(timecourse[TIME], timecourse[VALUE], label=well_name)

    # Label the axes and add a legend
    xlabel('Time') # TODO: automatically determine seconds or minutes, etc.
    ylabel('Value')
    fontP = FontProperties()
    fontP.set_size('small')
    legend(loc='upper center', prop=fontP, ncol=5, bbox_to_anchor=(0.5, 1.1),
         fancybox=True, shadow=True)

def plot_condition(wells, layout, condition_name):
    """Plot the replicates associated with a single experimental condition.

    Parameters
    ----------
    wells : dict
        The dict (returned by read_wallac) containing the well data.
    layout : dict specifying the plate layout.
        The keys are the names of the conditions; the values are lists
        of well names associated with that condition.
    condition_name : string
        The name of the condition to plot. Used to index into the plate
        layout dict.
    """
    wells_for_condition = layout[condition_name]
    for well_name in wells_for_condition:
        (time, value) = wells[well_name]
        plot(time, value)

def averages(wells, layout, stderr=False):
    # For every experimental condition...
    well_averages = collections.OrderedDict([])
    well_stds = collections.OrderedDict([])

    for condition_name, condition_list in layout.iteritems():
        num_conditions = len(condition_list)
        num_timepoints = len(wells[condition_list[0]][TIME])

        time_matrix = np.zeros((num_conditions, num_timepoints))
        value_matrix = np.zeros((num_conditions, num_timepoints))

        # Create a new entry in the averages dict of the form:
        #   well_averages[condition_name] = [[],[]]
        # Iterate over all of the wells
        for row_index, well_name in enumerate(condition_list):
            time_matrix[row_index,:] = wells[well_name][TIME]
            value_matrix[row_index,:] = wells[well_name][VALUE]

        # Calculate condition average over all wells
        time_averages = np.mean(time_matrix, axis=0)
        value_averages = np.mean(value_matrix, axis=0)
        time_stds = np.std(time_matrix, axis=0)
        if stderr:
            value_stds = np.std(value_matrix, axis=0) / np.sqrt(num_conditions)
        else:
            value_stds = np.std(value_matrix, axis=0)

        well_averages[condition_name] = [time_averages, value_averages]
        well_stds[condition_name] = [time_stds, value_stds]

    # end iteration over conditions
    return (well_averages, well_stds)

def distributions_by_condition(wells, layout):
    # For every experimental condition...
    distributions = collections.OrderedDict([])
    for condition_name, condition_list in layout.iteritems():
        num_conditions = len(condition_list)
        # Create a new entry in the averages dict consisting of a numpy
        # array
        distribution = []
        for row_index, well_name in enumerate(condition_list):
            distribution.append(wells[well_name])
        distributions[condition_name] = np.array(distribution)
    return distributions

def subtract_background(wells, background):
    """Subtract the given background vector (given as a single array)
    from each of the timecourse wells in the wells dict.
    """
    bgsub_timecourses = collections.OrderedDict([])
    for well_name in wells.keys():
        time = wells[well_name][TIME]
        tc = wells[well_name][VALUE]
        bgsub_tc = tc - background
        bgsub_timecourses[well_name] = [time, bgsub_tc]
    return bgsub_timecourses

def subtract_background_set(wells, background_dict):
    """Subtract the ordered dict of background vectors from the timecourse
    wells in the wells dict.
    """
    bgsub_timecourses = collections.OrderedDict([])
    for i, well_name in enumerate(wells.keys()):
        time = wells[well_name][TIME]
        tc = wells[well_name][VALUE]
        bgsub_tc = tc - background_dict.values()[i][VALUE]
        bgsub_timecourses[well_name] = [time, bgsub_tc]
    return bgsub_timecourses

def get_repeat_averages_by_well(wells):
    """For a set of replicate measurements made to establish a baseline
    (e.g., many measurements made of initial or final state)
    returns a dict mapping the well name to the value for that well,
    averaged over the serial replicates.
    """
    well_averages = collections.OrderedDict([])

    for well_name in wells.keys():
        repeat_measurements = wells[well_name][VALUE]
        well_averages[well_name] = np.mean(repeat_measurements)

    return well_averages

def get_first_points_by_well(wells):
    """From a timecourse, get the first timepoint (e.g., for use in
    normalization)."""
    first_points = collections.OrderedDict([])

    for well_name in wells.keys():
        repeat_measurements = wells[well_name][VALUE]
        first_points[well_name] = repeat_measurements[0]

    return first_points

def reset_first_timepoint_to_zero(wells):
    """From a timecourse, get the earliest timepoint and use that to subtract
    from all other timepoints."""
    min_time = np.inf
    for well_name in wells.keys():
        cur_min_time = np.min(wells[well_name][TIME])
        if cur_min_time < min_time:
            min_time = cur_min_time

    reset_wells = collections.OrderedDict([])
    for well_name in wells.keys():
        old_times = wells[well_name][TIME]
        new_times = [(old_time - min_time) for old_time in old_times]

        reset_wells[well_name] = [new_times, wells[well_name][VALUE]]

    return reset_wells

def reset_well_times(wells, offset=0):
    """For kinetics timecourses, resets the initial timepoint for each
    kinetics measurement to 0, or some offset."""
    reset_wells = collections.OrderedDict([])
    for well_name in wells.keys():
        min_time = np.min(wells[well_name][TIME])
        new_times = wells[well_name][TIME] - min_time + offset
        reset_wells[well_name] = [new_times, wells[well_name][VALUE]]

    return reset_wells

def get_baseline_value(baseline_wells, num_timepoints=10):
    """Takes a set of baseline wells and takes the first `num_timepoint`
    values, averages them for all of the baseline wells, and then returns
    the value, which can be used to normalize other wells for which a
    baseline measurement is missing or unreliable (e.g., in fast kinetics
    assays."""
    well_avgs = []
    for well_name in baseline_wells.keys():
        values = baseline_wells[well_name][VALUE]
        well_avgs.append(np.mean(values[0:num_timepoints]))
    well_avgs = np.array(well_avgs)
    return np.mean(well_avgs)

def get_normalized_well_timecourses(wells, initial_vals, final_vals):
    """For every well, normalizes to a percentage between the initial value and
    the final value."""
    normalized_wells = collections.OrderedDict([])

    for well_name in wells.keys():
        # If the initial vals are a dict, get the initial value from
        # the dict; otherwise, treat initial_vals as a fixed number offset
        try:
            initial_val = initial_vals[well_name]
        except IndexError:
            initial_val = initial_vals

        range = final_vals[well_name] - initial_val
        times = wells[well_name][TIME]
        values = wells[well_name][VALUE]
        normalized_values = (values - initial_val) / range
        normalized_wells[well_name] = [times, normalized_values]

    return normalized_wells

def get_average_pore_timecourses(normalized_wells):
    average_pores = collections.OrderedDict([])

    for well_name in normalized_wells.keys():
        times = normalized_wells[well_name][TIME]
        values = normalized_wells[well_name][VALUE]
        average_pore_values = -np.log(1 - values)
        average_pores[well_name] = [times, average_pore_values]

    return average_pores

def to_dataframe(mean_dict, sd_dict=None):
    """Convert timecourses to a pandas dataframe format.

    The keys in both dicts are expected to be space-separated strings, with the
    second item in the split string the Bax concentration.

    Parameters
    ----------
    mean_dict, sd_dict : dicts of lists
        dict mapping concentrations/conditions to the two-element list
        with a time vector and the fluorescence means or sds.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with the Bax concentration as the columns and a
        MultiIndexes on the rows with two levels: ("Timepoint", "Datatype").
        `Timepoint` is an integer index running from 0. `Datatype` is either
        'TIME', 'MEAN', or 'SD'. The entries are therefore either time values,
        mean fluorescence values, or the standard deviation of fluorescence
        values, depending on what type of row it is.
    """

    # Assume that all conditions have the same number of timepoints
    data_arr = None
    bax_concs = []
    row_tuples = None

    for conc_index, conc_string in enumerate(mean_dict.keys()):
        # Add this Bax concentration to the column index
        float_conc = float(conc_string.split(' ')[1])
        bax_concs.append(float_conc)

        # Get the values from the dicts
        time = mean_dict[conc_string][TIME]
        mean = mean_dict[conc_string][VALUE]
        if sd_dict is not None:
            sd = sd_dict[conc_string][VALUE]
        else:
            sd = 0

        if data_arr is None:
            # Initialize the size of the flattened data array
            data_arr = np.zeros([len(time) * 3, len(mean_dict.keys())])
        if row_tuples is None:
            # Initialize row_tuples for timepoint index
            row_tuples = []
            for i in range(len(time)):
                for row_tuple in zip([i]*3, ['TIME', 'MEAN', 'SD']):
                    row_tuples.append(row_tuple)
        # Fill in the data array
        data_arr[0::3, conc_index] = time
        data_arr[1::3, conc_index] = mean
        data_arr[2::3, conc_index] = sd

    bax_concs = np.array(bax_concs)

    # Build indices and dataframe
    row_index = pd.MultiIndex.from_tuples(row_tuples,
                                          names=('Timepoint', 'Datatype'))
    col_index = pd.Index(bax_concs, name='Bax')
    df = pd.DataFrame(data_arr, index=row_index, columns=col_index)

    return df

def truncate_timecourses(wells, start, end=None):
    trunc = collections.OrderedDict()

    for well in wells.keys():
        trunc[well] = []
        trunc[well].append(wells[well][TIME][start:end])
        trunc[well].append(wells[well][VALUE][start:end])
    return trunc

def get_replicates_for_condition(wells, layout, condition):
    """Takes a string representing a condition (e.g., 'Bax 100 nM'), gets the
    set of wells associated with the condition from the layout dict, then
    returns a time vector and a matrix with shape (num_timepoints,
    num_replicates) containing the timecourses for each replicate for
    that condition."""
    well_names = layout[condition]
    num_replicates = len(well_names)
    time_arr = wells[well_names[0]][TIME]
    num_timepoints = len(time_arr)
    replicate_matrix = np.zeros((num_timepoints, num_replicates))
    for i, well_name in enumerate(well_names):
        replicate_matrix[:,i] = wells[well_name][VALUE]
    return (time_arr, replicate_matrix)

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict

def wells_from_layout(layout):
    """Takes a layout dictionary and returns a list of the wells it contains.

    Examples
    --------

        >>> from collections import OrderedDict
        >>> layout = OrderedDict([('Bax 0 nM', ['A1', 'B1']),
        ...                       ('Bax 100 nM', ['A2', 'B2'])])
        >>> wells_from_layout(layout)
        ['A1', 'B1', 'A2', 'B2']
    """

    all_wells = []
    for label in layout.keys():
        cur_wells = layout[label]
        all_wells += cur_wells
    return all_wells

def dose_series_labels(reagent_name=None, initial_dose=None,
                       dilution_ratio=(2/3.),
                       num_doses=12, last_dose_zero=True,
                       units='nM', lowest_first=True):
    """Returns a list of strings describing dose-response conditions.

    Given a reagent name, and initial dose, a dilution ratio, and other
    parameters, returns a list of strings identifying the concentrations in
    each se of wells in the dilution series, e.g. ['Bax 400 nM', 'Bax 200 nM',
    'Bax 100 nM' ...] etc. This streamlines specifying the concentration
    conditions in different wells when plotting or analyzing data from
    plate-based assays.

    Parameters
    ----------
    reagent_name : string
        Name of the protein/compound in the dose series.
    initial_dose : float
        The highest dose in the dilution series.
    dilution_rate : float
        The scaling on the concentration going from the higher to lower
        concenrations in the series. So a 1:2 dilution series is 0.5;
        a 2:3 dilution series is 0.666.
    num_doses : int
        The number of doses in the dilution series, including the zero
        dose if present.
    last_dose_zero : boolean
        If True, the dilution series runs from the first to the (num_doses -
        1)th well, and the num_doses well is set to 0. If False, the dilution
        series span all num_doses wells.
    units : string
        The units identifier to be appended to the label, e.g. 'nM' or 'ng/mL'.
    lowest_first : boolean
        If True, the lowest concentration dose in the series is first in the
        list. If False, the highest concentration is first.

    Returns
    -------
    list of strings
        A list of labels identifying the concentrations of each dose in the
        series.

    Examples
    --------

    Getting a list of Bax concentrations:

        >>> dose_series_labels('Bax', 590., 1/2., num_doses=4)
        ['Bax 0.0 nM', 'Bax 147.5 nM', 'Bax 295.0 nM', 'Bax 590.0 nM']
    """

    if reagent_name is None or initial_dose is None:
        raise ValueError('Please specify values for the reagent_name '
                         'and initial_dose arguments!')
    dose_labels = []
    # Fill in the doses from highest concentration to lowest
    for i in range(num_doses):
        # If the last dose is zero, specify that here
        if last_dose_zero and i == num_doses - 1:
            dose = 0.
        # Otherwise, specify the diluted dose
        else:
            dose = initial_dose * (dilution_ratio ** i)
        # Format the label
        dose_label = '%s %.1f %s' % (reagent_name, dose, units)
        # Add it to the list
        dose_labels.append(dose_label)
    # Right now, the highest dose is the first in the list. If we want the
    # lowest (e.g., 0) dose first in the list (which is the default) we reverse
    # the list here:
    if lowest_first:
        dose_labels.reverse()

    return dose_labels

def replicate_well_list(start_row=None, end_row=None, start_col=None,
                        end_col=None, group_by_col=True):
    """Returns groups of well names for replicate wells.

    Expects that replicate wells occupy a rectangular block of contiguous wells
    in the plate.

    Parameters
    ----------
    start_row : string
        The first row in the replicate well area, e.g. 'A'.
    end_row : string
        The first row in the replicate well area, e.g. 'C'.
    start_col : int
        The first column in the replicate well area, e.g. 1.
    end_col : int
        The last column in the replicate well area, e.g. 12.
    group_by_col : boolean
        If True (default), wells with the same column index but different row
        indicies represent replicates of the same condition, with conditions
        varying from column to column. If False, replicates are grouped by row.

    Returns
    -------
    list of lists of strings
        Each list in the outer list represents a group of wells with
        the same condition, e.g. [['A1', 'B1', 'C1'], ['A2', 'B2', 'C2']]

    Examples
    --------
    >>> replicate_well_list('A', 'C', 1, 3)
    [['A1', 'B1', 'C1'], ['A2', 'B2', 'C2'], ['A3', 'B3', 'C3']]
    """

    if None in [start_row, end_row, start_col, end_col]:
        raise ValueError('Non-default values must be provided for row/col '
                         'arguments.')
    replicate_wells = []
    for row_char_num in range(ord(start_row), ord(end_row)+1):
        well_list = ['%s%d' % (chr(row_char_num), col)
                     for col in range(start_col, end_col+1)]
        replicate_wells.append(well_list)
    if group_by_col:
        replicate_wells = [list(tup) for tup in zip(*replicate_wells)]
    return replicate_wells

def dose_series_replicate_list(reagent_name=None, initial_dose=None,
                               dilution_ratio=(2/3.),
                               num_doses=12, last_dose_zero=True,
                               units='nM', lowest_first=True,
                               start_row=None, end_row=None, start_col=None,
                               end_col=None, group_by_col=True):
    """Returns strings/lists specifying a conditions in a plate.

    The strings/lists can be used to construct an OrderedDict specifying
    the plate layout.

    Expects that replicate wells occupy a rectangular block of contiguous wells
    in the plate. Note also that dilutions are expected to run from left to
    right in the plate, so the highest dose will be associated with the lowest
    row or column index. In this context the lowest_first argument merely
    specifies how the list is sorted, rather than the association of the labels
    to the wells.

    Parameters
    ----------
    reagent_name : string
        Name of the protein/compound in the dose series.
    initial_dose : float
        The highest dose in the dilution series.
    dilution_rate : float
        The scaling on the concentration going from the higher to lower
        concenrations in the series. So a 1:2 dilution series is 0.5;
        a 2:3 dilution series is 0.666.
    num_doses : int
        The number of doses in the dilution series, including the zero
        dose if present.
    last_dose_zero : boolean
        If True, the dilution series runs from the first to the (num_doses -
        1)th well, and the num_doses well is set to 0. If False, the dilution
        series span all num_doses wells.
    units : string
        The units identifier to be appended to the label, e.g. 'nM' or 'ng/mL'.
    lowest_first : boolean
        If True, the lowest concentration dose in the series is first in the
        list. If False, the highest concentration is first.
    start_row : string
        The first row in the replicate well area, e.g. 'A'.
    end_row : string
        The first row in the replicate well area, e.g. 'C'.
    start_col : int
        The first column in the replicate well area, e.g. 1.
    end_col : int
        The last column in the replicate well area, e.g. 12.
    group_by_col : boolean
        If True (default), wells with the same column index but different row
        indicies represent replicates of the same condition, with conditions
        varying from column to column. If False, replicates are grouped by row.

    Returns
    -------
    list of tuples of type (string, [list of strings])

    Examples
    --------
    >>> import collections
    >>> layout_dict = dose_series_replicate_list('Bax', 590, (2/3.), num_doses=3, start_row='A', end_row='C', start_col=1, end_col=3)
    >>> collections.OrderedDict(layout_dict)
    OrderedDict([('Bax 0.0 nM', ['A3', 'B3', 'C3']), ('Bax 393.3 nM', ['A2', 'B2', 'C2']), ('Bax 590.0 nM', ['A1', 'B1', 'C1'])])
    """

    if group_by_col and num_doses != end_col - start_col + 1:
        raise ValueError('The number of doses does not match the column range '
                         'in the plate layout.')

    rep_list = zip(dose_series_labels(reagent_name, initial_dose,
                                      dilution_ratio, num_doses, last_dose_zero,                                      units, lowest_first=False),
                   replicate_well_list(start_row, end_row, start_col, end_col,
                                       group_by_col))
    if lowest_first:
        rep_list.reverse()
    return rep_list

def plot_endpoints_vs_dose(wells, layout):
    endpoints = []
    conc_list = []
    for conc_str in layout.keys():
        conc = float(conc_str.split(' ')[1])
        conc_list.append(conc)
        endpoint = wells[conc_str][VALUE][-1]
        endpoints.append(endpoint)
    figure()
    plot(conc_list, endpoints, linestyle='', marker='o')
    xlabel('Dose')
    ylabel('RFU')

    return (conc_list, endpoints)

def row_wells(row_name, max_col=12):
    """Returns the list of wells in a row.

    Parameters
    ----------
    row_name : string
        Name of row, e.g. 'C'.
    max_col : int
        Maximum column index to use (inclusive). Default is 12.

    Returns
    -------
    list of strings containing the well names for the row.

    Examples
    --------
    >> row_wells('C', 12)
    ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12']
    """

    return ['%s%s' % (row_name, i) for i in range(1, max_col+1)]

def add_offset_vector(wells, offset_vector):
    """Add an offset vector to the time coordinates of each well.

    Parameters
    ----------
    wells : dict
        A dict containing well names as keys and [[time], [value]] pairs
        (list consisting of two numpy.arrays) as values.
    offset_vector : dict
        A dict containing well names as keys (matching the keys in the wells
        argument) and numbers as values, specifying the time offset for
        each well.

    Returns
    -------
    A copy of the wells dict with appropriately offset time coordinates.
    """
    wells_copy = copy(wells)
    for well_name in wells_copy.keys():
        well = wells_copy[well_name]
        well[TIME] += offset_vector[well_name]
    return wells_copy
