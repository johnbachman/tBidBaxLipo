"""A script to parse the Bid-568 + DiD liposome FRET data directly from the
Excel spreadsheet sent by Justin Kale into Python data files."""

from openpyxl import load_workbook
import numpy as np
import collections
from os.path import dirname, abspath, join
import sys
import tbidbaxlipo.data

# Load the data
data_path = dirname(sys.modules['tbidbaxlipo.data'].__file__)
data_file = abspath(join(data_path, '140429_Bid_membrane_FRET.xlsx'))
wb = load_workbook(data_file)

# Layout of the Excel spreadsheet
TIME_ROW_INDEX = 42
WELL_NAME_COL_INDEX = 0
FIRST_COL_INDEX = 1
FIRST_ROW_INDEX = 44
LAST_COL_INDEX = 154
LAST_ROW_INDEX = 140

# Load the first worksheet
sheet = wb.worksheets[0]

# Load the time coordinates
time = np.array([cell.value for cell in
                 sheet.rows[TIME_ROW_INDEX][FIRST_COL_INDEX:LAST_COL_INDEX]])

# FDA: WT cBid titration, DiD lipos, 10 nM Bid 568
fda = collections.OrderedDict([
         ('cBid 1000 nM', ['A1', 'B1', 'C1']),
         ('cBid 500 nM', ['A2', 'B2', 'C2']),
         ('cBid 250 nM', ['A3', 'B3', 'C3']),
         ('cBid 125 nM', ['A4', 'B4', 'C4']),
         ('cBid 62.5 nM', ['A5', 'B5', 'C5']),
         ('cBid 31.25 nM', ['A6', 'B6', 'C6']),
         ('cBid 15.6 nM', ['A7', 'B7', 'C7']),
         ('cBid 7.81 nM', ['A8', 'B8', 'C8']),
         ('cBid 3.91 nM', ['A9', 'B9', 'C9']),
         ('cBid 1.95 nM', ['A10', 'B10', 'C10']),
         ('cBid 0.98 nM', ['A11', 'B11', 'C11']),
         ('cBid 0 nM', ['A12', 'B12', 'C12']),
     ])

# FA: WT cBid titration, DiD lipos (no cBid-568)
fa = collections.OrderedDict([
         ('cBid 1000 nM', ['D1']),
         ('cBid 500 nM', ['D2']),
         ('cBid 250 nM', ['D3']),
         ('cBid 125 nM', ['D4']),
         ('cBid 62.5 nM', ['D5']),
         ('cBid 31.25 nM', ['D6']),
         ('cBid 15.6 nM', ['D7']),
         ('cBid 7.81 nM', ['D8']),
         ('cBid 3.91 nM', ['D9']),
         ('cBid 1.95 nM', ['D10']),
         ('cBid 0.98 nM', ['D11']),
         ('cBid 0 nM', ['D12'])
     ])
# FD: WT cBid titration, cBid-568, unlabeled lipos
fd = collections.OrderedDict([
         ('cBid 1000 nM', ['E1', 'F1', 'G1']),
         ('cBid 500 nM', ['E2', 'F2', 'G2']),
         ('cBid 250 nM', ['E3', 'F3', 'G3']),
         ('cBid 125 nM', ['E4', 'F4', 'G4']),
         ('cBid 62.5 nM', ['E5', 'F5', 'G5']),
         ('cBid 31.25 nM', ['E6', 'F6', 'G6']),
         ('cBid 15.6 nM', ['E7', 'F7', 'G7']),
         ('cBid 7.81 nM', ['E8', 'F8', 'G8']),
         ('cBid 3.91 nM', ['E9', 'F9', 'G9']),
         ('cBid 1.95 nM', ['E10', 'F10', 'G10']),
         ('cBid 0.98 nM', ['E11', 'F11', 'G11']),
         ('cBid 0 nM', ['E12', 'F12', 'G12']),
     ])
# BG: WT cBid titration, unlabeled lipos
bg = collections.OrderedDict([
         ('cBid 1000 nM', ['H1']),
         ('cBid 500 nM', ['H2']),
         ('cBid 250 nM', ['H3']),
         ('cBid 125 nM', ['H4']),
         ('cBid 62.5 nM', ['H5']),
         ('cBid 31.25 nM', ['H6']),
         ('cBid 15.6 nM', ['H7']),
         ('cBid 7.81 nM', ['H8']),
         ('cBid 3.91 nM', ['H9']),
         ('cBid 1.95 nM', ['H10']),
         ('cBid 0.98 nM', ['H11']),
         ('cBid 0 nM', ['H12'])
     ])

# The timecourses for each well will go in here
timecourse_wells = collections.OrderedDict()

# Iterate over the rows to get the individual trajectories
for row in sheet.rows[FIRST_ROW_INDEX:LAST_ROW_INDEX]:
    well_name = row[WELL_NAME_COL_INDEX].value
    timecourse = [cell.value for cell in row[FIRST_COL_INDEX:LAST_COL_INDEX]]
    timecourse_wells[well_name] = [time, timecourse]

