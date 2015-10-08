from tbidbaxlipo.data.parse_bid_bim_nbd_release import nbd_residues
import os

#nbd_residues = ['3', '5']
confs = range(2, 6)
activators = ['Bid', 'Bim']
reps = range(1, 4)
extensions = ['.fits', '.confs', '.conv', '.tri.low', '.tri.med', '.tri.hi']
dir_prefix = '/_static/pt_data1_plots/'

# Text for index file
index_text =  'Fits to Multiconformation Models\n'
index_text += '================================\n\n'
index_text += '.. toctree::\n'
index_text += '    :maxdepth: 2\n\n'

# Iterate over the mutants
for nbd_residue in nbd_residues:
    if nbd_residue == 'WT':
        continue
    text =  'NBD-%sC-Bax\n' % nbd_residue
    text += '===============\n\n'
    # Iterate over the activators
    for activator in activators:
        for rep in reps:
            text += '%s, Rep %d\n' % (activator, rep)
            text += '-----------------\n\n'
            # Iterate over the conformations
            for num_confs in confs:
                text += '%d conformations\n' % num_confs
                text += '~~~~~~~~~~~~~~~~~~~~\n\n'
                filename_base = os.path.join(dir_prefix,
                                 'pt_data1_%s_NBD_%s_r%d_%dconfs.mcmc' % \
                                 (activator, nbd_residue, rep, num_confs))
                # Iterate over the file extensions/images
                for extension in extensions:
                    png_file = filename_base + extension + '.png'
                    pdf_file = filename_base + extension + '.pdf'
                    text += '.. image:: %s\n\n' % png_file
                    text += ':download:`(pdf) <%s>`\n\n' % pdf_file
    # Write to file
    with open('NBD-%sC-Bax.rst' % nbd_residue, 'w') as f:
        f.write(text)

    index_text += '    NBD-%sC-Bax\n' % nbd_residue

# Write index file
with open('index.rst', 'w') as f:
    f.write(index_text)
