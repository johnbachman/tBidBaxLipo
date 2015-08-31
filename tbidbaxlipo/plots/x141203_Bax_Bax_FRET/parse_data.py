import tbidbaxlipo.data
import sys
from os.path import dirname, abspath, join
from openpyxl import load_workbook
import numpy as np
import pandas as pd

data_path = dirname(sys.modules['tbidbaxlipo.data'].__file__)
data_file = abspath(join(data_path,
           '2014-12-3 - Bax-Bax FRET with NBD and TbDPA RLS 126C and 54C.xlsx'))

nbd_residues = ['126', '54'] # This must match the spreadsheet
datatypes_time = ['Time', 'Release', 'FRET', 'NBD']
activator = 'cBid'
rep_index = 1

FIRST_ROW_INDEX = 2
LAST_ROW_INDEX = 128 # 54C ends on index 127, so will have a lot of NaNs

wb = load_workbook(data_file)

sheet = wb.worksheets[0]

col_index = 0
col_tuples = []
data = []

for nbd_ix, nbd_residue in enumerate(nbd_residues):
    # Iterate over the Time, Release, FRET, and NBD columns
    time_vector = None
    for dtype_ix, dtype in enumerate(datatypes_time):
        # Iterate and fill the array of column values
        value_vector = []
        data_col = sheet.columns[col_index]\
                                [FIRST_ROW_INDEX:LAST_ROW_INDEX]
        # Filter out empty cells
        for cell in data_col:
            if (cell.value == None):
                value_vector.append(np.nan)
            else:
                value_vector.append(cell.value)

        # If this is the time column, get and save the time vector
        if dtype_ix == 0:
            time_vector = value_vector
        # Otherwise, get the data and add it along with a time column
        else:
            assert time_vector is not None
            time_tuple = (activator, dtype, nbd_residue, rep_index, 'TIME')
            value_tuple = (activator, dtype, nbd_residue, rep_index, 'VALUE')
            col_tuples.append(time_tuple)
            col_tuples.append(value_tuple)
            data.append(time_vector)
            data.append(value_vector)
        # Increment the index of the current column
        col_index += 1

# Create the Pandas dataframe
data_matrix = np.array(data)
col_index = pd.MultiIndex.from_tuples(col_tuples,
                    names=('Activator', 'Datatype', 'NBD Site', 'Replicate',
                           'Column'))
df = pd.DataFrame(data_matrix.T,
                  index=range(LAST_ROW_INDEX - FIRST_ROW_INDEX),
                  columns=col_index)

