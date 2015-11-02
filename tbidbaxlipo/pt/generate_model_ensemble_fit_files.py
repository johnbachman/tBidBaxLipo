import yaml
from itertools import product, izip
from tbidbaxlipo.models import one_cpt
from tbidbaxlipo.models.nbd import multiconf
from copy import copy
import sys
import os

def model_product(model_dict):
    """A function to take the cross product of all feature implementations and
    return the set of dicts specifying one implementation for each feature."""
    return [dict(zip(model_dict, x))
            for x in product(*model_dict.values())]

def generate_files(args, basename):
    # Get the dict specifying the model ensemble feature set
    m = args['model']

    # Get the set of dicts specifying each individual implementation
    m_ensemble = model_product(m)
    dependencies_list = []

    for m in m_ensemble:
        # Multiconf model, supercedes any other model features
        if 'multiconf' in m or 'multiconf_nbd_fret' in m:
            multiconf_type = 'multiconf' if 'multiconf' in m \
                                         else 'multiconf_nbd_fret'
            # Is the NBD data normalized?
            norm_data = m['normalized_nbd_data']
            # Number of confs and whether the model is reversible
            try:
                num_confs = int(m[multiconf_type])
                reversible = False
            # If it's not an int, assume it's two-element list with a flag
            # specifying a reversible model, indicated by 'rev'
            except TypeError:
                num_confs = int(m[multiconf_type][0])
                rev = m[multiconf_type][1]
                if rev == 'rev':
                    reversible = True
                else:
                    raise Exception('Unknown multiconf model flag %s' % rev)
            # We don't need to build the model just to get the model name
            if reversible:
                model_name = '%dconfsrev' % num_confs
            else:
                model_name = '%dconfs' % num_confs
        # Mechanistic model
        else:
            # If Bid doesn't get to the membrane, activation by bound Bid will
            # never happen
            if 'bidtranslocation' in m and 'activation' in m and \
               m['bidtranslocation'] == 0 and m['activation'] == 2:
                continue
            # If activation is pseudo-first order (not Bid-dependent) don't use
            # models where we waste steps on translocating Bid
            if 'bidtranslocation' in m and 'activation' in m and \
               m['activation'] == 1 and m['bidtranslocation'] != 0:
                continue
            # If we're monitoring NBD as resulting from a dimer, but
            # dimerization doesn't happen, then we'll get nothing
            if 'nbd' in m and 'dimerization' in m and \
               (m['nbd'] == 2 or m['nbd'] == 3 or m['nbd'] == 4) and \
               m['dimerization'] == 0:
                continue
            # If we're monitoring NBD as resulting from a tBid/Bax complex,
            # then make sure that tBid/Bax complex can form
            if 'nbd' in m and 'activation' in m and \
                m['nbd'] == 4 and m['activation'] == 1:
                continue
            # If Bid FRET involves a complex of tBid with both mBax and iBax
            # (implementation 2) then make sure we have three-step activation,
            # which allows tBid to bind to iBax
            if 'bidfret' in m and 'activation' in m and \
               m['bidfret'] == 2 and not m['activation'] == 3:
                continue
            # It doesn't matter which builder is instantiated here because
            # the build_model_from_dict function is in the superclass, core.
            bd = one_cpt.Builder()
            bd.build_model_from_dict(m)
            model_name = bd.model.name

        # Build up the new yaml dict by copying all fitting parameters over...
        yaml_dict = copy(args)
        # ...and then filling in the parameters specifying this particular model
        yaml_dict['model'] = m
        model_filename = '%s_%s' % (basename, model_name)
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

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print("Please specify the name of the .yaml file with the model "
              "features to fit.")
        sys.exit(1)

    # Open the master .yaml file specifying the model ensemble
    ens_filename = sys.argv[1]
    with open(ens_filename) as ens_file:
        args = yaml.load(ens_file)

    # Get the basename for the created files
    basedir = os.path.dirname(ens_filename)
    filename = os.path.basename(ens_filename)
    basename = filename.split('.')[0]
    basename = os.path.join(basedir, basename)

    generate_files(args, basename)
