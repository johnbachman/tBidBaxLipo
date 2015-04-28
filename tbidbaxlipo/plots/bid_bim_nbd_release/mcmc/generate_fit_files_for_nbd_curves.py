import os
import sys
import yaml
from tbidbaxlipo.data.parse_bid_bim_nbd_release import nbd_residues

# Arguments shared across all fits
args = {
    'model': {
        'multiconf': [2, 3, 4],
        'normalized_nbd_data': [False]},
    'model_observable': ['NBD'],
    'global_initial_conditions': {
        'c0_scaling': 1.0, },
    'global_params': 'all',
    'local_params': [],
    'ntemps': 20,
    'highest_temp': -4,
    'nwalkers': 200,
    'nburnin': 20,
    'nsample': 20,
    'thin': 1,
}

basedir = sys.argv[1]
output_filename_pattern = 'pt_data1_%s_NBD_%s_r%s.fit.ensemble'
dependencies_list = []

# Iterate over the activators
for activator in ['Bid', 'Bim']:
    # Iterate over the NBD residues
    for nbd_residue in ['54', '126']:
        # Skip the wild type curves since there is no NBD trace
        if nbd_residue == 'WT':
            continue
        # Iterate over the replicates
        for rep_ix, rep_num in enumerate(range(1, 4)):
            # Initialize the data portion of the dict
            data_args = {
              'initial_condition_var': None,
              'module': 'tbidbaxlipo.plots.bid_bim_nbd_release.preprocess_data'}
            data_args['time_var'] = 'time_%s_NBD_%s_r%s' % \
                                    (activator, nbd_residue, rep_num)
            data_args['data_var'] = 'data_%s_NBD_%s_r%s' % \
                                    (activator, nbd_residue, rep_num)
            data_args['data_sigma_var'] = 'data_sigma_%s_NBD_%s' % \
                                          (activator, nbd_residue)
            args['data'] = data_args
            # Create the YAML file
            output_filename = output_filename_pattern % \
                              (activator, nbd_residue, rep_num)
            output_filename = os.path.join(basedir, output_filename)
            with open(output_filename, 'w') as f:
                f.write(yaml.dump(args, default_flow_style=False))
            dependencies_list.append(output_filename)

# Now write the file with the dependencies of the overall target on the
# list of .mcmc files
deps_filename = os.path.join(basedir, 'pt_data1.deps.txt')
target_name = 'pt_data1'
with open(deps_filename, 'w') as deps_file:
    ens_fit_file_list = [fname for fname in dependencies_list]
    #base_target = os.path.basename(basedir) # Strip off the directory info
    deps_file.write('%s: ' % target_name) # Strip off the directory info
    deps_file.write(' '.join(ens_fit_file_list))

