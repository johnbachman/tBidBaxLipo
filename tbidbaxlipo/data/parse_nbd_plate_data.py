"""
A script to parse the NBD fluorescence data from the Excel spreadsheet sent by
Justin Kale into Python data files.
"""

from openpyxl import load_workbook
from pylab import mean
from numpy import nan
import numpy as np

data_file = 'Compiled - Bid curves for NBD mutants.xlsx'
FIRST_COL_INDEX = 2
#FIRST_ROW_INDEX = 3
FIRST_ROW_INDEX = 2
LAST_COL_INDEX = 156

wb = load_workbook(data_file)
sheet = wb.worksheets[0]

data = []
# Get data column by column
for row in sheet.rows[FIRST_ROW_INDEX:len(sheet.columns)]:
    data_row = row[FIRST_COL_INDEX:LAST_COL_INDEX]
    row_arr = []
    # Iterate and fill the array of column values
    for cell in data_row:
        if (cell.value == None):
            row_arr.append(nan)
        else:
            row_arr.append(cell.value)
            #print cell.value
    data.append(row_arr)

time = data[0]
#nbd3c = data[1:4]
nbd62c = data[1:4]
#nbd120c = data[7:10]
#nbd122c = data[10:13]
#nbd126c = data[13:16]

#import ipdb; ipdb.set_trace()

output_filename = 'nbd_plate_data.py'
output_file = open(output_filename, 'w')
output_file.write('from numpy import mean, nan, array\n\n')

output_file.write('time = ' + np.array(data[0]).__repr__() + '\n\n')
output_file.write('nbd5c = ' + np.array(data[1:5]).__repr__() + '\n\n')
output_file.write('nbd15c = ' + np.array(data[7:11]).__repr__() + '\n\n')
output_file.write('nbd36c = ' + np.array(data[13:17]).__repr__() + '\n\n')
output_file.write('nbd40c = ' + np.array(data[19:23]).__repr__() + '\n\n')
output_file.write('nbd47c = ' + np.array(data[25:29]).__repr__() + '\n\n')
output_file.write('nbd54c = ' + np.array(data[31:35]).__repr__() + '\n\n')
output_file.write('nbd62c = ' + np.array(data[37:41]).__repr__() + '\n\n')
output_file.write('nbd68c = ' + np.array(data[43:47]).__repr__() + '\n\n')
output_file.write('nbd79c = ' + np.array(data[49:53]).__repr__() + '\n\n')
output_file.write('nbd120c = ' + np.array(data[55:59]).__repr__() + '\n\n')
output_file.write('nbd122c = ' + np.array(data[61:65]).__repr__() + '\n\n')
output_file.write('nbd126c = ' + np.array(data[67:71]).__repr__() + '\n\n')
output_file.write('nbd175c = ' + np.array(data[73:77]).__repr__() + '\n\n')
output_file.write('nbd179c = ' + np.array(data[79:83]).__repr__() + '\n\n')
output_file.write('nbd188c = ' + np.array(data[85:89]).__repr__() + '\n\n')

output_file.close()

