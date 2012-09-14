from openpyxl import load_workbook
from pylab import mean
from numpy import nan

#data_file = 'baxnbd.xlsx'
data_file = 'bax_nbd62c.xlsx'
FIRST_COL_INDEX = 0
#FIRST_ROW_INDEX = 3
FIRST_ROW_INDEX = 1

wb = load_workbook(data_file)
sheet = wb.worksheets[0]

data = []
# Get data column by column
for col in sheet.columns[FIRST_COL_INDEX:len(sheet.columns)]:
    column = col[FIRST_ROW_INDEX:len(col)]
    col_arr = []
    # Iterate and fill the array of column values
    for cell in column:
        if (cell.value == None):
            col_arr.append(nan)
        else:
            col_arr.append(cell.value)

    data.append(col_arr)

time = data[0]
#nbd3c = data[1:4]
nbd62c = data[1:4]
#nbd120c = data[7:10]
#nbd122c = data[10:13]
#nbd126c = data[13:16]


output_filename = 'nbd_data_62c.py'
output_file = open(output_filename, 'w')
output_file.write('from numpy import mean, nan\n\n')

output_file.write('time = ' + time.__repr__() + '\n\n')
#output_file.write('nbd3c = ' + nbd3c.__repr__() + '\n\n')
output_file.write('nbd62c = ' + nbd62c.__repr__() + '\n\n')
#output_file.write('nbd120c = ' + nbd120c.__repr__() + '\n\n')
#output_file.write('nbd122c = ' + nbd122c.__repr__() + '\n\n')
#output_file.write('nbd126c = ' + nbd126c.__repr__() + '\n\n')

#output_file.write('mean3c = mean(nbd3c, 0)\n\n')
#output_file.write('mean62c = mean(nbd62c, 0)\n\n')
#output_file.write('mean120c = mean(nbd120c, 0)\n\n')
#output_file.write('mean122c = mean(nbd122c, 0)\n\n')
#output_file.write('mean126c = mean(nbd126c, 0)\n\n')

output_file.close()


