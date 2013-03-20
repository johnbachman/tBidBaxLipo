import glob
from bayessb.report import Report
import bayessb.report.reporters
import tbidbaxlipo.reporters
from tbidbaxlipo.reporters import knowledge, topology, experiment
from tbidbaxlipo.nbd_mcmc_pysb import import_mcmc_groups
import sys
import pickle

#chain_file_list = glob.glob('/Volumes/data/computation/Bachman/'
#                            'c62_Baxbh3_nlkl/taird*.mcmc')
chain_file_list = glob.glob('/Volumes/data/computation/Bachman/'
                            'c62_Baxbh3_nlkl/taird_c62_iBax_50000_s0.mcmc')
#chain_file_list = glob.glob('ta*.mcmc')

chain_dict = import_mcmc_groups(chain_file_list)

rep = Report(chain_dict, [topology, bayessb.report.reporters,
                          knowledge, experiment])
rep.write_html_table_with_links('htmlreport.html')

with open('report.rep', 'w') as f:
    pickle.dump(rep, f)
