import os
import sys
import yaml
from tbidbaxlipo.data.parse_bid_bim_nbd_release import nbd_residues

# Arguments shared across all fits
args = {
    'model': {
        'multiconf': [[3, 'rev']], # 2, 3, 4
        'normalized_nbd_data': [False]},
    'model_observable': ['NBD'],
    'global_initial_conditions': {},
    'local_initial_condition': None,
    'global_params': 'all',
    'local_params': [],
    'ntemps': 50,
    'highest_temp': -6,
    'nwalkers': 400, # 400
    'nburnin': 5000,
    'nsample': 100,
    'thin': 1,
}

basedir = sys.argv[1]
output_target_pattern = 'pt_data1_%s_NBD_%s_r%s'
output_filename_pattern = output_target_pattern + '.fit.ensemble'
dependencies_list = []

# Iterate over the activators
for activator in ['Bid']:
    # Iterate over the NBD residues
    for nbd_residue in ['54']:
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
            data_args['nbd_ubound'] = 'data_%s_NBD_%s_r%s_ubound' % \
                                (activator, nbd_residue, rep_num)
            data_args['nbd_lbound'] = 'data_%s_NBD_%s_r%s_lbound' % \
                                (activator, nbd_residue, rep_num)
            data_args['nbd_f0'] = 'data_%s_NBD_%s_r%s_f0' % \
                                (activator, nbd_residue, rep_num)
            args['data'] = data_args
            # Create the YAML file
            output_target = output_target_pattern % \
                                (activator, nbd_residue, rep_num)
            output_filename = output_target + '.fit.ensemble'
            output_filename = os.path.join(basedir, output_filename)
            with open(output_filename, 'w') as f:
                f.write(yaml.dump(args, default_flow_style=False))
            dependencies_list.append(output_target)

# Now write the file with the dependencies of the overall target on the
# list of .mcmc files
deps_filename = os.path.join(basedir, 'pt_data1.deps.txt')
target_name = 'pt_data1'
with open(deps_filename, 'w') as deps_file:
    #base_target = os.path.basename(basedir) # Strip off the directory info
    # First, specify that the overall target depends on all the sub-targets
    deps_file.write('%s: ' % target_name) # Strip off the directory info
    deps_file.write(' '.join(dependencies_list) + '\n\n')
    # Next, for each dependency, specify that it depends on the corresponding
    # .deps.txt file, and include it
    for dep in dependencies_list:
        dep_filename = os.path.join(basedir, '%s.deps.txt' % dep)
        deps_file.write('%s: %s\n' % (dep, dep_filename))
        deps_file.write('-include %s\n\n' % dep_filename)
