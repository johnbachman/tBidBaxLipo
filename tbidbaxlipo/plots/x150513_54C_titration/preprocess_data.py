import os
import sys
import numpy as np
import tbidbaxlipo.data
from tbidbaxlipo.util.calculate_error_variance import calc_err_var
from openpyxl import load_workbook
from matplotlib import pyplot as plt

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
data_file = os.path.abspath(os.path.join(data_path,
                                        '2015-05-13 - 54C titration.xlsx'))


wb = load_workbook(data_file, data_only=True)
sheet = wb.worksheets[0]

FIRST_DATA_ROW = 6
LAST_DATA_ROW = 146

REL_TIME_COL_500 = 2
REL_DATA_COL_500 = 3
NBD_TIME_COL_500 = 6
NBD_DATA_COL_500 = 7

REL_TIME_COL_500 = 2
REL_DATA_COL_500 = 3
NBD_TIME_COL_500 = 6
NBD_DATA_COL_500 = 7

def get_col(col_index):
    return np.array([cell.value for cell in
                     sheet.columns[col_index][FIRST_DATA_ROW:LAST_DATA_ROW]])

data_dict = {}
data_dict['rt_500'] = get_col(2)
data_dict['ry_500'] = get_col(3)
data_dict['nt_500'] = get_col(6)
data_dict['ny_500'] = get_col(7)

data_dict['rt_250'] = get_col(10)
data_dict['ry_250'] = get_col(11)
data_dict['nt_250'] = get_col(14)
data_dict['ny_250'] = get_col(15)

data_dict['rt_125'] = get_col(18)
data_dict['ry_125'] = get_col(19)
data_dict['nt_125'] = get_col(22)
data_dict['ny_125'] = get_col(23)

data_dict['rt_62'] = get_col(26)
data_dict['ry_62'] = get_col(27)
data_dict['nt_62'] = get_col(30)
data_dict['ny_62'] = get_col(31)

def plot_data():
    plt.figure()
    plt.plot(data_dict['nt_500'], data_dict['ny_500'], label='Bplt 500 nM')
    plt.plot(data_dict['nt_250'], data_dict['ny_250'], label='Bplt 250 nM')
    plt.plot(data_dict['nt_125'], data_dict['ny_125'], label='Bplt 125 nM')
    plt.plot(data_dict['nt_62'], data_dict['ny_62'], label='Bplt 62.5 nM')
    plt.legend(loc='upper right')
    plt.xlabel('Time (sec)')
    plt.ylabel('NBD F/F0')

if __name__ == '__main__':
    plt.ion()
    plot_data()
