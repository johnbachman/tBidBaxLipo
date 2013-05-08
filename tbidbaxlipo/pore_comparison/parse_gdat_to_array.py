import glob
from pysb import bng
from tbidbaxlipo.pore_comparison.run_n_cpt import model
from tbidbaxlipo.models.n_cpt import get_dye_release
import pickle
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.util import color_iter

gdat_files = glob.glob('simdata/*.gdat')

xrecs = []
dr_all = []
num_sims = len(gdat_files)
for gdat_file in gdat_files:
    print "Loading %s" % gdat_file
    ssa_result = bng._parse_bng_outfile(gdat_file)
    xrecs.append(ssa_result)
    dr_all.append(get_dye_release(model, 'pores', ssa_result))

#xall = np.zeros((len(xrecs), len(xrecs[0]), len(xrecs[0].dtype)))
#for i, x in enumerate(xrecs):
#    xall[i,:,:] = xrecs[i].view(float).reshape(len(xrecs[i]), -1)
xall = np.array([x.tolist() for x in xrecs])
x_std = np.recarray(xrecs[0].shape, dtype=xrecs[0].dtype,
                    buf=np.std(xall, 0))
x_avg = np.recarray(xrecs[0].shape, dtype=xrecs[0].dtype,
                    buf=np.mean(xall, 0))

with open('x_avg.pck', 'w') as f:
    pickle.dump(x_avg, f)
with open('x_std.pck', 'w') as f:
    pickle.dump(x_std, f)
with open('dr.pck', 'w') as f:
    pickle.dump(dr_all, f)
