import glob
from bayessb.report import Report
import bayessb.report.reporters
import tbidbaxlipo.reporters
from tbidbaxlipo.reporters import knowledge, topology, experiment
from tbidbaxlipo.mcmc.nbd_mcmc import import_mcmc_groups
import sys
import pickle

#chain_file_list = glob.glob('/Volumes/data/computation/Bachman/'
#                            'multitemp/tardt_c62_iBax_50000_3_s*.mcmc')
chain_file_list = glob.glob('/files/ImStor/sorger/data/computation/Bachman/'
                            'pt_taid_1e6/taid*.mcmc')
#chain_file_list = glob.glob('ta*.mcmc')
#chain_file_list = glob.glob('/home/jab69/projects/tBidBaxLipo/tbidbaxlipo/mcmc/mpi1e6/ta*.mcmc')
print chain_file_list

chain_dict = import_mcmc_groups(chain_file_list)

rep = Report(chain_dict, [topology, bayessb.report.reporters,
                          knowledge, experiment])
rep.write_html_table_with_links('index.html')

with open('report.rep', 'w') as f:
    pickle.dump(rep, f)
