import pickle
import pkgutil

try:
    means = pickle.loads(pkgutil.get_data('tbidbaxlipo.simdata.sim_test',
                                          'means.pck'))
    stds = pickle.loads(pkgutil.get_data('tbidbaxlipo.simdata.sim_test',
                                          'stds.pck'))
except IOError:
    pass
