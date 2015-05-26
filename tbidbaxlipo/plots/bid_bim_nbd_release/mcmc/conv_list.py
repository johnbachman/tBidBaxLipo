import os
import numpy as np
import cPickle
import sys
from tbidbaxlipo.data.parse_bid_bim_nbd_release import nbd_residues
import subprocess

def save_last_position(mcmc_filename, pos_filename):
    # Get the sampler
    with open(mcmc_filename) as f:
        print("Loading %s" % mcmc_filename)
        (gf, sampler) = cPickle.load(f)
        # Get last position
        last_pos = sampler.chain[:,:,-1,:]
    with open(pos_filename, 'w') as f:
        # Get a default random meed
        np.random.seed(1)
        rs = np.random.get_state()
        print("Writing position to %s" % pos_filename)
        cPickle.dump((last_pos, rs), f)

# Iterate over the activators
for activator in ['Bid', 'Bim']:
    # Iterate over the NBD residues
    for nbd_residue in nbd_residues:
        # Iterate over the replicates
        for rep_ix, rep_num in enumerate(range(1, 4)):
            # Iterate over confs
            for num_confs in range(2, 6):
                filename_base = 'pt_data1_%s_NBD_%s_r%s_%sconfs' % \
                                (activator, nbd_residue, rep_num, num_confs)
                fit_filename = filename_base + '.fit'
                mcmc_filename = filename_base + '.mcmc'
                out_filename = filename_base + '.fit.out'
                err_filename = filename_base + '.fit.err'
                pos_filename = filename_base + '.fit.pos'
                posx_filename = pos_filename + 'x'
                if not os.path.isfile(mcmc_filename):
                    continue
                # Read the outputfile
                print("-------")
                converged = False
                try:
                    #print("Opening %s" % out_filename)
                    with open(out_filename) as outfile:
                        lines = outfile.readlines()
                        # Search the file backwards
                        for i in xrange(1, len(lines) + 1):
                            if lines[-i].startswith('-- Converged!') or \
                               lines[-i].startswith('Passes:  True'):
                                converged = True
                                break
                        print("%s %s" %
                              (out_filename,
                               "Converged" if converged else "Not converged"))
                except:
                    print("Could not open %s" % out_filename)
                # Check for a position file
                if os.path.isfile(pos_filename):
                    print("Found %s" % pos_filename)
                    position_filename = pos_filename
                elif os.path.isfile(posx_filename):
                    print("Found %s" % posx_filename)
                    position_filename = posx_filename
                else:
                    print("Missing %s and %s" % (pos_filename, posx_filename))
                    save_last_position(mcmc_filename, posx_filename)
                # If not converged, submit
                if not converged:
                    if sys.argv[1] == 'qsub':
                        cmd_list = ['qsub', '-b', 'y', '-cwd', '-V', '-o',
                                    out_filename, '-e', err_filename, '-pe',
                                    'orte', '32', 'mpirun', 'python',
                                    '../../../pt/run_pt.py', fit_filename,
                                    '1', position_filename]
                        print ' '.join(cmd_list)
                        subprocess.call(cmd_list)
                    elif sys.argv[1] == 'bsub':
                        cmd_list = ['bsub', '-q', 'parallel', '-n', '10',
                                    '-W', '100:00', '-N', '-o',
                                    out_filename, '-a', 'openmpi', 'mpirun.lsf',
                                    'python', '-m', 'tbidbaxlipo.pt.run_pt',
                                    fit_filename, '1', position_filename]
                        print ' '.join(cmd_list)
                        subprocess.call(cmd_list)
                    else:
                        raise Exception("Invalid scheduler command")


