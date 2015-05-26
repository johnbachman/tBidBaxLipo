import os
import numpy as np
import cPickle
import sys

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
    for nbd_residue in ['36', '47', '54', '62', '68', '79', '126']: # 54, 126
        # Iterate over the replicates
        for rep_ix, rep_num in enumerate(range(1, 4)):
            # Iterate over confs
            for num_confs in range(2, 6):
                filename_base = 'pt_data1_%s_NBD_%s_r%s_%sconfs' % \
                                (activator, nbd_residue, rep_num, num_confs)
                mcmc_filename = filename_base + '.mcmc'
                out_filename = filename_base + '.fit.out'
                pos_filename = filename_base + '.fit.pos'
                posx_filename = pos_filename + 'x'
                # Read the outputfile
                try:
                    #print("Opening %s" % out_filename)
                    with open(out_filename) as outfile:
                        lines = outfile.readlines()
                        # Search the file backwards
                        converged = False
                        for i in xrange(1, len(lines) + 1):
                            if lines[-i].startswith('-- Converged!') or \
                               lines[-i].startswith('Passes:  True'):
                                converged = True
                                break
                        #print("%s %s" %
                        #      (out_filename,
                        #       "Converged" if converged else "Not converged"))
                except:
                    print("Could not open %s" % out_filename)
                # Check for a position file
                if not (os.path.isfile(pos_filename) or \
                        os.path.isfile(posx_filename)):
                    print("Missing %s and %s" % (pos_filename, posx_filename))
                    save_last_position(mcmc_filename, pos_filename + 'x')
                elif os.path.isfile(pos_filename):
                    print("Found %s" % pos_filename)
                elif os.path.isfile(posx_filename):
                    print("Found %s" % posx_filename)
                else:
                    assert False


