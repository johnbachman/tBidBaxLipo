from openpyxl import load_workbook
import pandas as pd
from itertools import product
import numpy as np
import tbidbaxlipo.data
import numpy as np

import os
import sys

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)

if True:
    data_file = os.path.abspath(os.path.join(data_path,
            '2017-09-19 - BAX-BAX FRET dataset V2 Compiled.xlsx'))
    # Zero-indexed
    FIRST_ROW_INDEX = 8
    LAST_ROW_INDEX = 134
else:
    data_file = os.path.abspath(os.path.join(data_path,
            '2016 - Bax-Bax FRET dataset.xlsx'))
    # Zero-indexed
    FIRST_ROW_INDEX = 3
    LAST_ROW_INDEX = 128

# These lists must match the order in the spreadsheet exactly
nbd_residues = ['3', '36', '54', '68', '120', '126', '138',
                '175', '179']
activators = ['Bid']
#datatypes = ['Release', 'FRET', 'NBD']
datatypes_time = ['Time', 'Release', 'FRET', 'NBD']
reps = [1, 2, 3]
col_types = ['TIME', 'VALUE']


wb = load_workbook(data_file, data_only=True)

# Load the first worksheet
sheet = wb.worksheets[0]

# Initialize some values before iterating
col_index = 0
col_tuples = []
data = []

activator = 'Bid'

# For all residues...
for nbd_ix, nbd_residue in enumerate(nbd_residues):
    # ...and replicates...
    for rep_ix, rep in enumerate(reps):
        time_vector = None
        # Now iterate over the Time, Release, FRET, NBD columns
        for dtype_ix, dtype in enumerate(datatypes_time):
            # Iterate and fill the array of column values
            vector = []
            data_col = sheet.columns[col_index]\
                                     [FIRST_ROW_INDEX:LAST_ROW_INDEX]
            for cell in data_col:
                if (cell.value == None):
                    vector.append(np.nan)
                else:
                    vector.append(float(cell.value))
            if dtype_ix == 0:
                time_vector = vector
            else:
                assert time_vector is not None
                time_tuple = (activator, dtype, nbd_residue, rep, 'TIME')
                value_tuple = (activator, dtype, nbd_residue, rep, 'VALUE')
                col_tuples.append(time_tuple)
                col_tuples.append(value_tuple)
                data.append(time_vector)
                data.append(vector)
            col_index += 1
            # Increment the index of the current column

# Create the Pandas dataframe
data_matrix = np.array(data)
col_multi_index = pd.MultiIndex.from_tuples(col_tuples,
                    names=('Activator', 'Datatype', 'NBD Site', 'Replicate',
                           'Column'))

df = pd.DataFrame(data_matrix.T,
                  index=range(LAST_ROW_INDEX - FIRST_ROW_INDEX),
                  columns=col_multi_index)
