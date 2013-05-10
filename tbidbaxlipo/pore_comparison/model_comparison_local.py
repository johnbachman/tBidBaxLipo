from pylab import figure, title, legend
from tbidbaxlipo.models import one_cpt, site_cpt, n_cpt
from tbidbaxlipo.pore_comparison.run_n_cpt import params_dict, tmax, model
from tbidbaxlipo.models import n_cpt
import pickle

"""
params_dict = {'Vesicles_0':1, 'tBid_0':1, 'Bax_0':50,
'Bax_transloc_kr':1e-2}

params_dict = {
        'Vesicles_0': 1,
        'tBid_0': 1,
        'Bax_0': 50,
        'tBid_transloc_kf': 0.1,
        'tBid_transloc_kr': 0.1,
        'Bax_transloc_kf': 0.01,
        'Bax_transloc_kr': 0.01,
        'tBid_mBax_kf': 0.01,
        'tBid_mBax_kr': 1.5,
        'tBid_iBax_kc': 0.1
}
"""
params_dict = {
    #'c62_scaling': 1.4022073763353842,
    'Vesicles_0': 50,
    'tBid_0': 0,
    'Bax_0': 50,
    #'tBid_transloc_kf': 0.1,
    #'tBid_transloc_kr': 0.87463035606442208,
    #'tBid_transloc_kr': 0,
    'Bax_transloc_kf': 0.01,
    'Bax_transloc_kr': 0.1,
    #'tBid_mBax_kf': 0.0010752348795729722,
    #'tBid_mBax_kr': 0.00021614462019141619,
    #'tBid_iBax_kc': 0.031624502179168776,
    'pore_formation_rate_k': 1e-3,
    #'tBid_iBax_kf': 0.02455155520493121,
    #'tBid_iBax_kr': 0.0023148402339379995,
    #'Bax_dimerization_kf': 0.00087114734334575424,
    #'Bax_dimerization_kr': 0.0012875415792703037,
}

sc = 2

# Compartment-based vs. ODE
b1c = one_cpt.Builder(params_dict=params_dict)
b1c.build_model_bax_schwarz()
b1c.run_model(tmax=tmax)

#x_avg = pickle.load(open('x_avg.pck'))
#x_std = pickle.load(open('x_std.pck'))
#dr = pickle.load(open('dr.pck'))
#num_sims = len(dr)
#print "Num sims %d" % num_sims
#n_cpt.plot_simulation(x_avg, x_std, dr, model, num_sims)

bnc = n_cpt.Builder(scaling_factor=sc, params_dict=params_dict)
bnc.build_model_bax_schwarz()
bnc.run_model(tmax=tmax, num_sims=3)

title('ODE vs. Compartment-based')

