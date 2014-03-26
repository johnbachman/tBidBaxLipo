from openpyxl import load_workbook
import pandas as pd
from itertools import product
import numpy as np
import tbidbaxlipo.data

import os
import sys

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
data_file = os.path.abspath(os.path.join(data_path,
            'Compiled Release percentages and NBD F-F0.xlsx'))

nbd_residues = ['WT', '3', '5', '15', '36', '40', '47', '54', '62', '68', '79',
                '120', '122', '126', '138', '151', '175', '179', '184', '188']
activators = ['Bid', 'Bim']
datatypes = ['Release', 'NBD']
reps = [1, 2, 3]
col_types = ['TIME', 'VALUE']

col_tuples = list(product(activators, datatypes, nbd_residues, reps, col_types))

# Zero-indexed
TIME_COL_INDEX = 0
FIRST_COL_INDEX = 1
FIRST_ROW_INDEX = 1
LAST_COL_INDEX = 61
LAST_ROW_INDEX = 137

wb = load_workbook(data_file)

# Load the first worksheet
sheet0 = wb.worksheets[0]


def get_data_from_sheet(sheet, activator, datatype):
    # Get data from time col
    # Get data column by column as a list of lists
    num_reps_per_mutant = 3
    time_vector = [cell.value for cell in
            sheet.columns[0][FIRST_ROW_INDEX:LAST_ROW_INDEX]]

    data = []
    tuples = []
    for i, col in enumerate(sheet.columns[FIRST_COL_INDEX:]):
        # The replicate index
        rep_index = (i % num_reps_per_mutant) + 1
        # Get the name for this mutant--it's in the first col for the rep,
        # in the first row
        if rep_index == 1:
            mutant_name = str(col[0].value)

        # Make the tuples that we'll need to create the dataframe
        time_tuple = (activator, datatype, mutant_name, rep_index, 'TIME')
        value_tuple = (activator, datatype, mutant_name, rep_index, 'VALUE')

        # First, add the time vector for the data on this sheet
        tuples.append(time_tuple)
        data.append(time_vector)

        data_col = col[FIRST_ROW_INDEX:LAST_ROW_INDEX]
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

(bid_release_tuples, bid_release_data) = \
                    get_data_from_sheet(wb.worksheets[0], 'Bid', 'Release')
(bim_release_tuples, bim_release_data) = \
                    get_data_from_sheet(wb.worksheets[1], 'Bim', 'Release')
(bid_nbd_tuples, bid_nbd_data) = \
                    get_data_from_sheet(wb.worksheets[2], 'Bid', 'NBD')
(bim_nbd_tuples, bim_nbd_data) = \
                    get_data_from_sheet(wb.worksheets[3], 'Bim', 'NBD')

col_tuples = bid_release_tuples + bim_release_tuples + \
             bid_nbd_tuples + bim_nbd_tuples
data_matrix = np.concatenate((bid_release_data, bim_release_data,
                             bid_nbd_data, bim_nbd_data))

col_index = pd.MultiIndex.from_tuples(col_tuples,
                    names=('Activator', 'Datatype', 'NBD Site', 'Replicate', 'Column'))

df = pd.DataFrame(data_matrix.T,
                  index=range(LAST_ROW_INDEX - FIRST_ROW_INDEX),
                  columns=col_index)

