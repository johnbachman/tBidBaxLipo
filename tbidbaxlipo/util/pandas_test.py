import pandas as pd
import numpy as np

lipo_concs = [2, 4]
bax_concs = [0, 5, 10]

num_timepoints = 3

# Build the row index
row_tuples = []
for j in range(num_timepoints):
    row_tuples.append((j, 'TIME'))
    row_tuples.append((j, 'MEAN'))
    row_tuples.append((j, 'SD'))

# Build the column index
col_tuples = []
for lipo_conc in lipo_concs:
    for bax_conc in bax_concs:
        col_tuples.append((lipo_conc, bax_conc))

# Build the dataset
num_columns = len(col_tuples)
data = np.zeros([3*num_timepoints, num_columns])
for i in range(num_columns): # Fill in the columns (concentrations)
    for j in range(num_timepoints): # Fill in the rows (timepoints)
        data[(j*3)+0, i] = j+1
        data[(j*3)+1, i] = np.log(j+1)
        data[(j*3)+2, i] = 0.5

# Build the dataframe
row_index = pd.MultiIndex.from_tuples(row_tuples, names=('Timepoint', 'Datatype'))
col_index = pd.MultiIndex.from_tuples(col_tuples, names=('Liposomes', 'Bax'))
df = pd.DataFrame(data, index=row_index, columns=col_index)
print df

# Get timecourse for a given concentration
tc = df[(4, 5)]

# Get the first timepoint:
tc[0]

# Now, I want to get all of the fluorescence values across the timecourse.
# Interestingly, these two approaches don't work, since 'MEAN'
# is from the "inner" level; slicing is not automatic/intelligent:
# tc['MEAN']
# tc.ix['MEAN']

# However, these approaches work:
print tc[:, 'MEAN']
print tc.ix[:, 'MEAN']
