"""
A set of functions useful for analyzing and plotting data from timecourse
assays from the Wallac plate reader.

- Export the data, sorted by well. To do this, take the Microsoft Excel output
  file, go to the "List" tab; highlight the data; then go to the "Data" menu
  and select sort by "Well". Finally export the data to CSV as a "Windows
  Comma Separated (.csv)" file.

  NOTE: the csv reader module used in this script appears to have a problem
  with .csv files exported from Excel using the default .csv settings. However,
  it works fine with files exported using the .csv "Windows" setting.

- Create a plate layout dict, for example, the following:

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

"""
# TODO Write a function that will read the initial and final .csv or
# excel files and use the values to normalize the timecourses

import csv
import datetime as dt
#from pylab import *
from matplotlib.pyplot import *
from matplotlib.font_manager import FontProperties

# Indices of the exported CSV
WELL_COL = 2
TIME_COL = 4
VALUE_COL = 5

# We store the well data in a list of lists; the first list
# has the time coordinates, the second list has the values
# (see read_wallac, below)
TIME = 0
VALUE = 1

def read_wallac(csv_file):
    """Parses the plate reader data from a CSV file.
  
    NOTE: The data must be sorted by well before exporting to CSV. See the
    instructions in the module docstring for more details.
  
    Parameters
    ----------
    csv_file : string
        The name of the CSV file to open.

    Returns
    -------
    A dict mapping the name of the well to a list of lists. For each dict
    entry, the first 
    element in the outer list is a list of the time coordinates associated
    with each timepoint; the second element in the outer list is the list of
    assay values for each timepoint.

    For example, wells['A01'] returns

        [[79.8, 145.43, ..., 2639.35, 2704.98],
         [2285,   2318, ...,    2323,    2304]]
    """

    # The plan: iterate over the CSV file, collecting time coordinates
    # and values until we find that we're at a new well; at which point
    # we update the value of cur_well, make a new entry in the dict,
    # and continue.
    csv_reader = csv.reader(open(csv_file, 'r'))
  
    # Initialize the empty dict
    wells = {}
    # We're not currently on any well!
    cur_well = ""
    # Iterate over the rows:
    for row in csv_reader:
        well_name = row[WELL_COL]

        # If we're at the header line, skip it
        if well_name == 'Well':
           continue
        # If it's a new well, update cur_well and initialize the dict with
        # a pair of empty lists
        elif cur_well != well_name:
            cur_well = well_name
            wells[well_name] = [np.array([]),np.array([])]
        # If it's a well we've seen before, add the time/val pair to the dict
        else:
            time = dt.datetime.strptime(row[TIME_COL], '%H:%M:%S.%f')
            secs = (time.hour * 3600) + (time.minute * 60) + time.second + \
                   (time.microsecond * 0.000001)
            value = int(row[VALUE_COL])
            wells[cur_well][TIME] = np.append(wells[cur_well][TIME], secs)
            wells[cur_well][VALUE] = np.append(wells[cur_well][VALUE], value)

    return wells

def plot_all(wells, errors=None):
    """Plots all of the timecourses in the dict.

    Parameters
    ----------
    wells : dict of timecourse data of the type returned by read_wallac.
    """

    for wellname, wellval in iter(sorted(wells.iteritems())):
        if errors is None:
            plot(wellval[TIME], wellval[VALUE], label=wellname)
        else:
            errorbar(wellval[TIME], wellval[VALUE],
                    yerr=errors[wellname][VALUE], label=wellname)

    # Label the axes and add a legend
    xlabel('Time') # TODO: automatically determine seconds or minutes, etc.
    ylabel('Value')
    fontP = FontProperties()
    fontP.set_size('small')
    legend(loc='upper center', prop=fontP, ncol=5, bbox_to_anchor=(0.5, 1.1),
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

def plot_condition(wells, layoujt, condition_name):
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

"""
def averages(wells, layout):
    # Calculates the average, at each timepoint, across replicates for each
    # experimental condition described in the plate layout.
    #

    # For every experimental condition...
    well_averages = {}
    for condition_name, condition_list in layout.iteritems():
        # Create a new entry in the averages dict of the form:
        #   well_averages[condition_name] = [[],[]]
        # Iterate over all of the wells
        counter = 0
        for well_name in condition_list:
            if counter == 0:
                time_sum = wells[well_name][TIME]
                value_sum = wells[well_name][VALUE]
            else:
                time_sum = map(lambda x, y: x + y, time_sum,
                               wells[well_name][TIME])
                value_sum = map(lambda x, y: x + y, value_sum,
                                wells[well_name][VALUE])
            counter += 1

        # Calculate condition average over all wells
        time_average = map(lambda x: float(x) / float(counter), time_sum)
        value_average = map(lambda x: float(x) / float(counter), value_sum)
        well_averages[condition_name] = [time_average, value_average]

    # end iteration over conditions
    return well_averages
"""

def averages(wells, layout):
    # For every experimental condition...
    well_averages = {}
    well_stds = {}
    
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
        value_stds = np.std(value_matrix, axis=0)

        well_averages[condition_name] = [time_averages, value_averages]
        well_stds[condition_name] = [time_stds, value_stds]

    # end iteration over conditions
    return (well_averages, well_stds)

def get_repeat_averages_by_well(wells):
    """For a set of replicate measurements made to establish a baseline
    (e.g., many measurements made of initial or final state)
    returns a dict mapping the well name to the value for that well,
    averaged over the serial replicates.
    """
    well_averages = {}

    for well_name in wells.keys():
        repeat_measurements = wells[well_name][VALUE]
        well_averages[well_name] = np.mean(repeat_measurements)

    return well_averages

def get_normalized_well_timecourses(wells, initial_vals, final_vals):
    """For every well, normalizes to a percentage between the initial value and
    the final value."""
    normalized_wells = {}

    for well_name in wells.keys():
        range = final_vals[well_name] - initial_vals[well_name]
        times = wells[well_name][TIME]
        values = wells[well_name][VALUE]
        normalized_values = (values - initial_vals[well_name]) / range
        normalized_wells[well_name] = [times, normalized_values]
    
    return normalized_wells

