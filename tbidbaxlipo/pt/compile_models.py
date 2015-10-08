"""Load all of the .yaml files in this directory, create models from them,
and give a pilot run to pre-compile the RHS functions."""

from pysb.integrate import Solver
import numpy as np
import yaml
import sys
from tbidbaxlipo.models import one_cpt
from tbidbaxlipo.models.nbd import multiconf
from pysb.integrate import Solver

Solver._use_inline = True

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
        # Number of confs and whether the model is reversible
        try:
            num_confs = int(args['model']['multiconf'])
            reversible = False
        # If it's not an int, assume it's two-element list with a flag
        # specifying a reversible model, indicated by 'rev'
        except TypeError:
            num_confs = int(args['model']['multiconf'][0])
            rev = args['model']['multiconf'][1]
            if rev == 'rev':
                reversible = True
            else:
                raise Exception('Unknown multiconf model flag %s' % rev)

        norm_data = args['model']['normalized_nbd_data'] = True
        bd.build_model_multiconf(num_confs, 1, normalized_data=norm_data,
                                 reversible=reversible)
    else:
        bd = one_cpt.Builder()
        bd.build_model_from_dict(args['model'])

    t = np.linspace(0, 1e4, 1e3)
    sol = Solver(bd.model, t)
    sol.run()
