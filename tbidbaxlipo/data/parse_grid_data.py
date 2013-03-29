"""
A script to parse the dye release data from the Excel spreadsheets sent by
David Andrews into Python data files.
"""

import os
import glob
import re
from openpyxl import load_workbook

dir_name = '../../../Data From David Andrews/Gridv2/Graphs'

FIRST_COL_INDEX = 1 # Note zero indexing!
HEADER_ROW_INDEX = 1 # Note zero indexing!
data_range = 'B2:N16'

time = []
grid = {}
mt1_grid = {}

for data_file in glob.glob(os.path.join(dir_name, 'Graphs*.xlsx')):
    print "Parsing file " + data_file + "..."

    # Parse out the amount of liposomes from the filename
    # Note that data_file contains the full (absolute) path, not just the filename
    m = re.search('Graphs_([0-9.]*) uL_Lipos.xlsx', data_file)
    #m = re.search('Graphs_([0-9.]*) uL Lipos.xlsx', data_file)
    if(m):
        lipo_conc = m.group(1)
        grid[lipo_conc] = {}
        mt1_grid[lipo_conc] = {}
    else:
        raise Exception("Couldn't parse liposome concentration.")

    # Initialize no cBid (Bax only) sub-array to an empty dict
    grid[lipo_conc]['0'] = {}
    mt1_grid[lipo_conc]['0'] = {}

    # Initialize liposome only sub-array (for this liposome concentration) to an empty list
    grid[lipo_conc]['0']['0'] = []
    mt1_grid[lipo_conc]['0']['0'] = []

    # Iterate over worksheets
    wb = load_workbook(data_file)
    for sheet in wb.worksheets:

        # Parse the worksheet title to see which one we're in
        bax_match = re.search('Bax', sheet.title)
        wt_cBid_match = re.search('([0-9.]*) nM wt', sheet.title)
        mt1_cBid_match = re.search('([0-9.]*) nM mt1', sheet.title)

        # Skip the worksheets that are for fixed concentrations of Bax,
        # since they contain both the wt cBid and the cBid mt dose-responses;
        # easier to parse if we stick to the cBid worksheets
        if (bax_match): 
            continue
        # Also skip the mt1 cBid titrations, at least for now
        elif (mt1_cBid_match):
            continue
        # Parse out the cBid concentration           
        elif (wt_cBid_match):
            cBid_conc = wt_cBid_match.group(1)
            grid[lipo_conc][cBid_conc] = {}

            # Get the cells column by column
            data_columns = sheet.columns[FIRST_COL_INDEX:len(sheet.columns)]
            for col in data_columns:
                timecourse = col[HEADER_ROW_INDEX:len(col)]
                header_string = ''
                arr_to_fill = None

                for cell in timecourse:
                    # A bit wacky, but works to identify the header row
                    if (cell.row == HEADER_ROW_INDEX + 1): 
                        header_string = cell.value

                        ## Now that we've got the header, make sure that we put the data
                        ## values in the right place:

                        # Check and see if this is the time column 
                        if (re.match('Time \(min\)', header_string)):
                            if (not time):
                                arr_to_fill = time
                            else:
                                break
                        # Check and see if this is the lipo only column
                        elif (re.match('Liposomes', header_string)):
                            if (not grid[lipo_conc]['0']['0']):
                                arr_to_fill = grid[lipo_conc]['0']['0']
                            else:
                                break
                        # Check and see if this is a Bax only (no cBid) column
                        elif (not re.search('cBid', header_string)):
                            # Parse out the Bax concentration
                            m = re.match('([0-9.]*) nM Bax', cell.value)
                            if (not m):
                                raise Exception('Error parsing the Bax concentration for cell ' + str(cell.value))
                            else:
                                bax_conc = m.group(1)
                                grid[lipo_conc]['0'][bax_conc] = []
                                arr_to_fill = grid[lipo_conc]['0'][bax_conc]
                        # Otherwise, this should be a cBid column, so we should parse it
                        else:
                            # Parse out the Bax concentration
                            # Two conditions to check for here--those with Bax and those without (cBid only)
                            # -- first, the cBid only condition
                            if (re.search('cBid', cell.value) and not re.search('Bax', cell.value)):
                                grid[lipo_conc][cBid_conc]['0'] = []
                                arr_to_fill = grid[lipo_conc][cBid_conc]['0']
                            # -- second, the case with cBid and Bax
                            else:
                                # Parse out the Bax concentration
                                m = re.search('cBid \+ ([0-9.]*) nM Bax', cell.value)
                                if (m):
                                    bax_conc = m.group(1)
                                    grid[lipo_conc][cBid_conc][bax_conc] = []
                                    arr_to_fill = grid[lipo_conc][cBid_conc][bax_conc]
                                else:
                                    raise Exception('Error parsing the Bax concentration for cell ' + str(cell.value))
                        # end header parsing
                    # end header row condition
                    ## If we're not in the header row, the target for the data
                    ## values should have been set, so just stick the value in
                    else:
                        arr_to_fill.append(cell.value)
                # end timecourse iteration
            # end column (concentration) iteration
        # end wt cBid worksheet parsing
        else:
            print('Worksheet title did not contain Bax, wt cBid, or mt1 cBid!')
    # end worksheet (cBid conc) iteration
# end file (lipo conc) iteration

output_filename = 'grid_data_v2.py'
output_file = open(output_filename, 'w')
output_file.write('time_min = ' + time.__repr__() + '\n\n')
output_file.write('time = [min*60 for min in time_min]\n\n')
output_file.write('data_tbidmaj = ' + grid.__repr__())
output_file.close()

