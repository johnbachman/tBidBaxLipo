import os
import numpy as np
import cPickle
import sys
from tbidbaxlipo.data.parse_bid_bim_nbd_release import nbd_residues
import subprocess

usage_msg =  "Usage:\n"
usage_msg += " python submit_chain_plots.py job_scheduler output_dir file_glob"

if len(sys.argv) < 3:
    print(usage_msg)
    sys.exit(1)

job_scheduler = sys.argv[1]
output_dir = sys.argv[2]

# Iterate over the activators
files = sys.argv[3:]
for mcmc_filename in files:
    out_filename = mcmc_filename + '.plots.out'
    err_filename = mcmc_filename + '.plots.err'

    if job_scheduler == 'qsub':
        cmd_list = ['qsub', '-b', 'y', '-cwd', '-V', '-o',
                    out_filename, '-e', err_filename,
                    'python', '-m', 'tbidbaxlipo.pt.show_chain',
                    mcmc_filename, output_dir]
        print ' '.join(cmd_list)
        subprocess.call(cmd_list)
    elif job_scheduler == 'bsub':
        cmd_list = ['bsub', '-q', 'short',
                    '-W', '1:00', '-N', '-o',
                    out_filename,
                    'python', '-m', 'tbidbaxlipo.pt.show_chain',
                    mcmc_filename, output_dir]
        print ' '.join(cmd_list)
        subprocess.call(cmd_list)
    else:
        raise Exception("Invalid scheduler command")


