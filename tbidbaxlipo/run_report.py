import glob
from bayessb.report import Report
import bayessb.report.reporters
import tbidbaxlipo.reporters
from tbidbaxlipo.nbd_mcmc_pysb import import_mcmc_groups
import sys

print sys.argv[1]
chain_file_list = glob.glob(sys.argv[1])
#chain_file_list = glob.glob('/files/ImStor/sorger/data/computation/Bachman/' +
#                            'tard_c62_*.mcmc')
chain_dict = import_mcmc_groups(chain_file_list)

rep = Report(chain_dict, [bayessb.report.reporters, tbidbaxlipo.reporters])
rep.write_html_table_with_links('htmlreport.html')


