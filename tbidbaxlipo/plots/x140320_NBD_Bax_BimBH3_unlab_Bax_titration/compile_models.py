"""Load all of the .yaml files in this directory, create models from them,
and give a pilot run to pre-compile the RHS functions."""

from pysb.integrate import Solver
import numpy as np
import yaml
import glob

from tbidbaxlipo.models.one_cpt import Builder

yaml_filenames = glob.glob('./*.yaml')

for yaml_filename in yaml_filenames:
    with open(yaml_filename) as yaml_file:
        args = yaml.load(yaml_file)

    bd = Builder()
    bd.build_model_from_dict(args['model'])
    t = np.linspace(0, 1e4, 1e3)
    sol = Solver(bd.model, t)
    sol.run()
