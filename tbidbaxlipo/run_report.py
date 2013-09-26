"""
A script for generating reports on MCMC fits of tBid/Bax/Lipo models
to data.

Usage::

    python -m tbidbaxlipo.run_report *.mcmc

The glob specifies which pickled MCMC files are to be included. To change the
reports that are run, edit this script.
"""

import glob
from bayessb.report import Report
import bayessb.report.reporters
from tbidbaxlipo.reporters import knowledge, topology, experiment, residuals, \
                                  titration_fits
from tbidbaxlipo.mcmc import import_mcmc_groups
import sys
import pickle

# Get the list of MCMC files
if len(sys.argv) <= 1:
    print("The list of MCMC files (e.g. as a glob) must be included.")
    sys.exit()
chain_file_list = sys.argv[1:]
print chain_file_list

# Parse the list of MCMC files into a dict of lists
chain_dict = import_mcmc_groups(chain_file_list)

# Run the reporters
rep = Report(chain_dict, [bayessb.report.reporters, titration_fits])
#rep = Report(chain_dict, [bayessb.report.reporters, residuals])
#rep = Report(chain_dict, [topology, bayessb.report.reporters,
#                          knowledge, experiment])

# Write the HTML table with supporting files
rep.write_html_table_with_links('index.html')

# Pickle the report object
with open('report.rep', 'w') as f:
    pickle.dump(rep, f)

