import glob
from bayessb.report import Report
import bayessb.report.reporters
import tbidbaxlipo.reporters
from tbidbaxlipo.reporters import knowledge, topology
from tbidbaxlipo.nbd_mcmc_pysb import import_mcmc_groups
import sys

chain_file_list = glob.glob('/Volumes/data/computation/Bachman/c3_testrun/' +
                            'tard_c3_*.mcmc')
chain_file_list = glob.glob('/Volumes/data/computation/Bachman/c3_testrun/' +
                            'tard_c3_iBax_50000_s9.mcmc')
#chain_file_list = glob.glob('/files/ImStor/sorger/data/computation/Bachman/' +
#                            'tard_c62_*.mcmc')
chain_dict = import_mcmc_groups(chain_file_list)

#rep = Report(chain_dict, [bayessb.report.reporters])
rep = Report(chain_dict, [bayessb.report.reporters,
                         topology, knowledge])
rep.write_html_table_with_links('htmlreport.html')

