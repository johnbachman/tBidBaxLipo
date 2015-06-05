from tbidbaxlipo.models.one_cpt import Builder
from matplotlib import pyplot as plt
from preprocess_data import time_36
from pysb.integrate import Solver

# Note activation 3 does not go with Bax reversal 1
# (unless reversal is rewritten to be independent of binding state)

model_dict = {
    'bidtranslocation': 1,
    'baxtranslocation': 1,
    'activation': 3,
    'reversal': 0,
    'nbd': 5,
    'bidfret': 1,
    'baxfret': 0,
    'bleach': 0,
    }

params_dict = {
        'Vesicles_0': 1.55,
        'Bax_0': 0.,
        'Bax_NBD_0': 100.,
        'tBid_Bax_ins_bh3_kf': 0.01,
        'tBid_Bax_ins_bh3_kr': 0.15,
        #'tBid_Bax_ins_bh3_kc': 0.01,
        #'tBid_Bax_mem_bh3_kr': 0.01,
}

bd = Builder(params_dict=params_dict)
bd.build_model_from_dict(model_dict)

t = time_36
s = Solver(bd.model, t)
s.run()

plt.ion()
plt.figure('NBD')
plt.plot(t, s.yexpr['NBD'])
plt.figure('FRET')
plt.plot(t, s.yexpr['BidFRET'])


