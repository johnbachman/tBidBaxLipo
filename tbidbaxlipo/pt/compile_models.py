"""Load all of the .yaml files in this directory, create models from them,
and give a pilot run to pre-compile the RHS functions."""

from pysb.integrate import Solver
import numpy as np
import yaml
import sys
from tbidbaxlipo.models import one_cpt
from tbidbaxlipo.models.nbd import multiconf

if len(sys.argv) == 1:
    print("Please specify one or model .yaml files with model information "
          "for pre-compilation.")
    sys.exit(1)

yaml_filenames = sys.argv[1:]

for yaml_filename in yaml_filenames:
    with open(yaml_filename) as yaml_file:
        args = yaml.load(yaml_file)
    # Check for whether this is a multiconf model or not
    if 'multiconf' in args['model']:
        bd = multiconf.Builder()
        num_confs = args['model']['multiconf']
        norm_data = args['model']['normalized_nbd_data'] = True
        bd.build_model_multiconf(num_confs, 1, normalized_data=norm_data)
    else:
        bd = one_cpt.Builder()
        bd.build_model_from_dict(args['model'])

    t = np.linspace(0, 1e4, 1e3)
    sol = Solver(bd.model, t)
    sol.run()
