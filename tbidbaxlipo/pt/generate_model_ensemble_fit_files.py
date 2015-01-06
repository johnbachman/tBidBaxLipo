import yaml
from itertools import product, izip
from tbidbaxlipo.models.one_cpt import Builder
from copy import copy
import sys
import os

if len(sys.argv) == 1:
    print("Please specify the name of the .yaml file with the model "
          "features to fit.")
    sys.exit(1)

# Open the master .yaml file specifying the model ensemble
ens_filename = sys.argv[1]
with open(ens_filename) as ens_file:
    args = yaml.load(ens_file)

# Get the dict specifying the model ensemble feature set
m = args['model']

# A function to take the cross product of all feature implementations and
# return the set of dicts specifying one implementation for each feature
def model_product(model_dict):
    return (dict(izip(model_dict, x))
            for x in product(*model_dict.itervalues()))

# Get the set of dict specifying each individual implementation
m_ensemble = model_product(m)
model_names = []
dependencies_list = []
basename = ens_filename.split('.')[0]
for m in m_ensemble:
    # Build the model from the dict
    bd = Builder()
    bd.build_model_from_dict(m)
    # Build up the new yaml dict by copying all fitting parameters over...
    yaml_dict = copy(args)
    # ...and then filling in the parameters specifying this particular model
    yaml_dict['model'] = m
    model_filename = '%s_%s' % (basename, bd.model_name)
    with open('%s.fit' % model_filename, 'w') as output_file:
        output_file.write(yaml.dump(yaml_dict, default_flow_style=False))
    dependencies_list.append(model_filename)

# Now write the file with the dependencies of the overall target on the
# list of .mcmc files
deps_filename = '%s.deps.txt' % basename
with open(deps_filename, 'w') as deps_file:
    mcmc_file_list = ['%s.mcmc' % fname for fname in dependencies_list]
    base_target = os.path.basename(basename) # Strip off the directory info
    deps_file.write('%s: ' % base_target) # Strip off the directory info
    deps_file.write(' '.join(mcmc_file_list))

