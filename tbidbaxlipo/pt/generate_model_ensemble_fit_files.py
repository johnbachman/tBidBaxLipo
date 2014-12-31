import yaml
from itertools import product, izip
from tbidbaxlipo.models.one_cpt import Builder
from copy import copy
import sys

if len(sys.argv) == 1:
    print("Please specify the name of the .yaml file with the model "
          "features to fit.")
    sys.exit(1)

yaml_filename = sys.argv[1]
with open(yaml_filename) as yaml_file:
    args = yaml.load(yaml_file)

m = args['model']

def model_product(model_dict):
    return (dict(izip(model_dict, x))
            for x in product(*model_dict.itervalues()))

m_ensemble = model_product(m)
model_names = []
for m in m_ensemble:
    bd = Builder()
    bd.build_model_from_dict(m)
    yaml_dict = copy(args)
    yaml_dict['model'] = m
    basename = yaml_filename.split('.')[0]
    filename = '%s_%s.yaml' % (basename, bd.model_name)
    with open(filename, 'w') as output_file:
        output_file.write(yaml.dump(yaml_dict, default_flow_style=False))
