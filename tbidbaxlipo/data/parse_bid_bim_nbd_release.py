from openpyxl import load_workbook
import pandas as pd
from itertools import product
import numpy as np
import tbidbaxlipo.data

import os
import sys

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
data_file = os.path.abspath(os.path.join(data_path,
              '2015-03-06-Compiled Release percentages and NBD F-F0_V2.xlsx'))
            #'Compiled Release percentages and NBD F-F0.xlsx'))

nbd_residues = ['WT', '3', '5', '15', '36', '40', '47', '54', '62', '68', '79',
                '120', '122', '126', '138', '151', '175', '179', '184', '188']
activators = ['Bid', 'Bim']
datatypes = ['Release', 'NBD']
reps = [1, 2, 3]
col_types = ['TIME', 'VALUE']

# Zero-indexed
TIME_COL_INDEX = 0
FIRST_COL_INDEX = 1
LAST_COL_INDEX = 61
NUM_ROWS = 135 # Actually 136, but the Bim sheet is missing a row

def get_data_from_sheet(sheet, activator, datatype, first_row_index, num_rows,
                        mutant_name_row):
    last_row_index = first_row_index + num_rows
    # Get data from time col
    # Get data column by column as a list of lists
    num_reps_per_mutant = 3
    time_vector = [cell.value for cell in
            sheet.columns[0][first_row_index:last_row_index]]

    data = []
    tuples = []
    for i, col in enumerate(sheet.columns[FIRST_COL_INDEX:]):
        # The replicate index
        rep_index = (i % num_reps_per_mutant) + 1
        # Get the name for this mutant--it's in the first col for the rep,
        # in the first row
        if rep_index == 1:
            mutant_name = str(col[mutant_name_row].value)

        # Make the tuples that we'll need to create the dataframe
        time_tuple = (activator, datatype, mutant_name, rep_index, 'TIME')
        value_tuple = (activator, datatype, mutant_name, rep_index, 'VALUE')

        # First, add the time vector for the data on this sheet
        tuples.append(time_tuple)
        data.append(time_vector)

        data_col = col[first_row_index:last_row_index]
        value_vector = []
        # Iterate and fill the array of column values
        for cell in data_col:
            if (cell.value == None):
                value_vector.append(nan)
            else:
                value_vector.append(cell.value)
        tuples.append(value_tuple)
        data.append(value_vector)
    return [tuples, np.array(data)]

def get_labeling_ratios(sheet, mutant_name_row, ratio_row):
    ratios = {}
    for i, col in enumerate(sheet.columns[FIRST_COL_INDEX:]):
        mutant_name = str(col[mutant_name_row].value)
        if mutant_name == '':
            continue
        labeling_ratio = col[ratio_row].value
        if mutant_name not in ratios:
            ratios[mutant_name] = labeling_ratio
    return ratios

def load_data():
    global df, labeling_ratios

    col_tuples = list(product(activators, datatypes, nbd_residues, reps,
                              col_types))


    wb = load_workbook(data_file, data_only=True)

    (bid_release_tuples, bid_release_data) = \
                        get_data_from_sheet(wb.worksheets[0], 'Bid', 'Release',
                                            1, NUM_ROWS, 0)
    (bim_release_tuples, bim_release_data) = \
                        get_data_from_sheet(wb.worksheets[1], 'Bim', 'Release',
                                            1, NUM_ROWS, 0)

    # Get data for the raw (uncorrected) values
    (bid_nbd_tuples, bid_nbd_data) = \
                        get_data_from_sheet(wb.worksheets[4], 'Bid', 'NBD',
                                            3, NUM_ROWS, 2)
    (bim_nbd_tuples, bim_nbd_data) = \
                        get_data_from_sheet(wb.worksheets[5], 'Bim', 'NBD',
                                            3, NUM_ROWS, 2)

    # Get labeling ratios
    labeling_ratios = get_labeling_ratios(wb.worksheets[4], 2, 1)

    col_tuples = bid_release_tuples + bim_release_tuples + \
                 bid_nbd_tuples + bim_nbd_tuples
    data_matrix = np.concatenate((bid_release_data, bim_release_data,
                                 bid_nbd_data, bim_nbd_data))

    col_index = pd.MultiIndex.from_tuples(col_tuples,
                        names=('Activator', 'Datatype', 'NBD Site', 'Replicate',
                               'Column'))

    df = pd.DataFrame(data_matrix.T, index=range(NUM_ROWS), columns=col_index)

# Load the data
load_data()
