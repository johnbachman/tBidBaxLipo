"""
A script to parse the NBD fluorescence data from the Excel spreadsheet sent by
Justin Kale into Python data files.
"""

from openpyxl import load_workbook
import numpy as np
import pandas as pd
from itertools import product
import pickle

data_file = '09-26-2012-Bax Liposome titration kinetics.xlsx'
# Zero-indexed
FIRST_COL_INDEX = 1
FIRST_ROW_INDEX = 2
LAST_COL_INDEX = 13
LAST_ROW_INDEX = 362

wb = load_workbook(data_file)

# Load the first worksheet
sheet0 = wb.worksheets[0]

# Load the time coordinates
time = np.array([cell.value for cell in
                 sheet0.columns[0][FIRST_ROW_INDEX:LAST_ROW_INDEX]])
# Convert to seconds
time = time * 60

# Load the liposome concentrations
lipo_concs = np.array([cell.value for cell in
                      sheet0.rows[1][FIRST_COL_INDEX:LAST_COL_INDEX]])

def get_data_from_sheet(sheet):
    # Get data column by column as a list of lists
    data = []
    for row in sheet.rows[FIRST_ROW_INDEX:LAST_ROW_INDEX]:
        data_row = row[FIRST_COL_INDEX:LAST_COL_INDEX]
        row_arr = []
        # Iterate and fill the array of column values
        for cell in data_row:
            if (cell.value == None):
                row_arr.append(nan)
            else:
                row_arr.append(cell.value)
        data.append(row_arr)
    return np.array(data)

# Get the data from each sheet
bax50nm = get_data_from_sheet(sheet0)
bax100nm = get_data_from_sheet(wb.worksheets[1])
bax200nm = get_data_from_sheet(wb.worksheets[2])
data_matrix = np.concatenate((bax50nm, bax100nm, bax200nm), axis=1)

# Build the column index
col_tuples = list(product([50., 100., 200.], lipo_concs))
col_index = pd.MultiIndex.from_tuples(col_tuples, names=('Bax', 'Liposomes'))

# Build the dataset
data = pd.DataFrame(data_matrix, index=time, columns=col_index)

if __name__ == '__main__':
    # Pickle it
    pickle.dump(data, open('nbd_c126_titration_data.pck', 'w'))
