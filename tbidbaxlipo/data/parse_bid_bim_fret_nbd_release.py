from openpyxl import load_workbook
import pandas as pd
from itertools import product
import numpy as np
import tbidbaxlipo.data
import numpy as np

import os
import sys

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
data_file = os.path.abspath(os.path.join(data_path,
        '2015-02-27 - Bid-Bim FRET with Bax NBD mutants in Tb-DPA Liposomes compiled reps+WT Bax.xlsx'))

# These lists must match the order in the spreadsheet exactly
nbd_residues = ['WT', '3', '15', '36', '47', '54', '62', '122', '126', '138',
                '151', '175', '184']
activators = ['Bid', 'Bim']
datatypes_time = ['Time', 'Release', 'FRET', 'NBD']
reps = [1, 2, 3]
col_types = ['TIME', 'VALUE']

# Zero-indexed
FIRST_ROW_INDEX = 4
LAST_ROW_INDEX = 128

wb = load_workbook(data_file)

# Load the first worksheet
sheet = wb.worksheets[0]

# Initialize some values before iterating
col_index = 0
col_tuples = []
data = []

# For all residues...
for nbd_ix, nbd_residue in enumerate(nbd_residues):
    # ...and replicates...
    for rep_ix, rep in enumerate(reps):
        # ...and activators...
        for act_ix, activator in enumerate(activators):
            # Set the time vector to None so we can make sure that it is set
            # correctly
            time_vector = None
            # Now iterate over the Time, Release, FRET, NBD columns
            for dtype_ix, dtype in enumerate(datatypes_time):
                # Skip the FRET and NBD datatypes if we're currently on the
                # WT Bax section
                if nbd_residue == 'WT' and (dtype == 'FRET' or dtype == 'NBD'):
                    col_index += 1
                    continue

                # Iterate and fill the array of column values
                value_vector = []
                data_col = sheet.columns[col_index]\
                                         [FIRST_ROW_INDEX:LAST_ROW_INDEX]
                for cell in data_col:
                    if (cell.value == None):
                        value_vector.append(np.nan)
                    else:
                        value_vector.append(float(cell.value))

                # If this is the time column, get and save the time vector
                if dtype_ix == 0:
                    time_vector = value_vector
                # Otherwise, get the data and add it along with a time column
                else:
                    assert time_vector is not None
                    time_tuple = (activator, dtype, nbd_residue, rep, 'TIME')
                    value_tuple = (activator, dtype, nbd_residue, rep, 'VALUE')
                    col_tuples.append(time_tuple)
                    col_tuples.append(value_tuple)
                    data.append(time_vector)
                    data.append(value_vector)
                # Increment the index of the current column
                col_index += 1

# Create the Pandas dataframe
data_matrix = np.array(data)
col_multi_index = pd.MultiIndex.from_tuples(col_tuples,
                    names=('Activator', 'Datatype', 'NBD Site', 'Replicate',
                           'Column'))
#import ipdb; ipdb.set_trace()
df = pd.DataFrame(data_matrix.T,
                  index=range(LAST_ROW_INDEX - FIRST_ROW_INDEX),
                  columns=col_multi_index)
