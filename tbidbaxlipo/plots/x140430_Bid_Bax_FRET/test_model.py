from tbidbaxlipo.models.one_cpt import Builder
from matplotlib import pyplot as plt
from preprocess_data import time_36
from pysb.integrate import Solver
import numpy as np

# Note activation 3 does not go with Bax reversal 1
# (unless reversal is rewritten to be independent of binding state)

model_dict = {
    'bidtranslocation': 1,
    'baxtranslocation': 1,
    'activation': 2,
    'reversal': 0,
    'nbd': 5,
    'bidfret': 1,
    'baxfret': 0,
    'bleach': 0,
    }

#params = [-1.89187047, -4.22923156, -5.47509186, -4.50648375,
#          -5.07143377,  3.70925668, -1.46660821, 0.982916,
#          0.32449226, 1.92229914]

def get_params_dict(params):

    names = [
         'tBid_transloc_kf', 'tBid_transloc_kr', 'Bax_transloc_kf',
        'Bax_transloc_kr', 'tBid_Bax_mem_bh3_kf', 'tBid_Bax_mem_bh3_kr',
        'tBid_Bax_ins_bh3_kc', 'c1_scaling', 'c2_scaling', 'bid_fret1',
    ]
    tuples = zip(names, params)
    params_dict = dict(tuples)
    params_dict['Vesicles_0'] = 1.55
    params_dict['Bax_0'] = 0.
    params_dict['Bax_NBD_0'] = 100.
    params_dict['tBid_0'] = 20.
    return params_dict

params_dict= dict([('Vesicles_0', 1.55),
 ('tBid_0', 20.0),
 ('Bax_0', 0.0),
 ('Bax_DAC_0', 0.0),
 ('Bax_NBD_0', 100.0),
 ('tBid_transloc_kf', 0.00030071703231223111),
 ('tBid_transloc_kr', 0.00015915898888252953),
 ('Bax_transloc_kf', 5.2689491196789301e-05),
 ('Bax_transloc_kr', 0.037348626553698416),
 ('tBid_Bax_mem_bh3_kf', 14.441045395282226),
 ('tBid_Bax_mem_bh3_kr', 0.0035791420490433657),
 ('tBid_Bax_ins_bh3_kc', 0.012599318998715796),
 ('c0_scaling', 1.0),
 ('c1_scaling', 2.1056228028869444),
 ('c2_scaling', 1.594115169777732),
 ('bid_fret1', 60.736662816869426),])

"""
params_dict = dict([
 ('Vesicles_0', 1.55),
 ('tBid_0', 20.0),
 ('Bax_0', 0.0),
 ('Bax_DAC_0', 0.0),
 ('Bax_NBD_0', 100.0),
 ('tBid_transloc_kf', 0.004750801171922365),
 ('tBid_transloc_kr', 0.34500310052186267),
 ('Bax_transloc_kf', 0.00045733459025325604),
 ('Bax_transloc_kr', 0.16429023567939888),
 ('tBid_Bax_mem_bh3_kf', 0.25107314636248984),
 ('tBid_Bax_mem_bh3_kr', 104.74686479435127),
 ('tBid_Bax_ins_bh3_kc', 0.10737077147536374),
 ('c0_scaling', 1.0),
 ('c1_scaling', 5.4631866488434095),
 ('c2_scaling', 6.2661077850658842),
 ('bid_fret1', 43.855089658242584), ])
"""

"""
params_dict = {
        'Vesicles_0': 1.55,
        'Bax_0': 0.,
        'Bax_NBD_0': 100.,
        'tBid_0': 20.,
        'tBid_transloc_kf': 10**-1.89187047,
        'tBid_transloc_kr': 10**-4.22923156,
        'Bax_transloc_kf': 10**-5.47509186,
        'Bax_transloc_kr': 10**-4.50648375,
        'tBid_Bax_mem_bh3_kf': 10**-5.07143377,
        'tBid_Bax_mem_bh3_kr': 10**3.70925668,
        'tBid_Bax_ins_bh3_kc': 10**-1.46660821,
        'c1_scaling': 10**0.982916,
        'c2_scaling': 10**0.32449226,
        'bid_fret1': 10**1.92229914,
}
"""

bd = Builder(params_dict=params_dict)
bd.build_model_from_dict(model_dict)
print(bd.model.parameters)
t = time_36
s = Solver(bd.model, t)
s.run()

n_species = s.y.shape[1]

plt.ion()
for spec_ix in range(n_species):
    plt.figure()
    plt.title(bd.model.species[spec_ix])
    plt.plot(t, s.y[:,spec_ix])

#plt.ion()
#plt.figure('NBD')
#plt.plot(t, s.yexpr['NBD'])
#plt.figure('FRET')
#plt.plot(t, s.yexpr['BidFRET'])


