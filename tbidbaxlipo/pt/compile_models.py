"""Load all of the .yaml files in this directory, create models from them,
and give a pilot run to pre-compile the RHS functions."""

from pysb.integrate import Solver
import numpy as np
import yaml
import sys
from tbidbaxlipo.models.one_cpt import Builder

if len(sys.argv) == 1:
    print("Please specify one or model .yaml files with model information "
          "for pre-compilation.")
    sys.exit(1)

yaml_filenames = sys.argv[1:]

for yaml_filename in yaml_filenames:
    with open(yaml_filename) as yaml_file:
        args = yaml.load(yaml_file)

    bd = Builder()
    bd.build_model_from_dict(args['model'])
    t = np.linspace(0, 1e4, 1e3)
    sol = Solver(bd.model, t)
    sol.run()
