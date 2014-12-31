import yaml
from itertools import product, izip
from tbidbaxlipo.models.one_cpt import Builder

yaml_filename = 'model_range.yaml'
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
    model_names.append(bd.model_name)
