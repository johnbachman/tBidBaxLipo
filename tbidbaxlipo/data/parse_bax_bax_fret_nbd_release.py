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
        '2016 - Bax-Bax FRET dataset.xlsx'))

# These lists must match the order in the spreadsheet exactly
nbd_residues = ['3', '36', '54', '68', '126', '134', '120',
                '175', '179']
activators = ['Bid']
datatypes = ['Release', 'FRET', 'NBD']
reps = [1, 2, 3]
col_types = ['TIME', 'VALUE']

# Zero-indexed
FIRST_ROW_INDEX = 3
# Should be 128, but chopped off the last few points to avoid massive outliers
LAST_ROW_INDEX = 128

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
        # Now iterate over the Release, FRET, NBD columns
        for dtype_ix, dtype in enumerate(datatypes):
            for col_type_ix, col_type in enumerate(col_types):
                # Iterate and fill the array of column values
                vector = []
                data_col = sheet.columns[col_index]\
                                         [FIRST_ROW_INDEX:LAST_ROW_INDEX]
                factor = 100 if ((dtype == 'Release' or dtype == 'FRET') \
                                 and col_type == 'VALUE') else 1
                for cell in data_col:
                    if (cell.value == None):
                        vector.append(np.nan)
                    else:
                        vector.append(float(cell.value) * factor)
                col_tuples.append((activator, dtype, nbd_residue, rep,
                                   col_type))
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
