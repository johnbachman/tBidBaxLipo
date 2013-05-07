from pylab import figure, title, legend
from tbidbaxlipo.models import one_cpt, site_cpt, n_cpt

tmax = 500

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
    'Vesicles_0': 5,
    'tBid_0': 5,
    'Bax_0': 100,
    'tBid_transloc_kf': 0.1,
    'tBid_transloc_kr': 0.87463035606442208,
    #'tBid_transloc_kr': 0,
    'Bax_transloc_kf': 0.01,
    'Bax_transloc_kr': 1.8033442952454604e-05,
    'tBid_mBax_kf': 0.0010752348795729722,
    'tBid_mBax_kr': 0.00021614462019141619,
    'tBid_iBax_kc': 0.031624502179168776,
    'pore_formation_rate_k': 1e-3,
    #'tBid_iBax_kf': 0.02455155520493121,
    #'tBid_iBax_kr': 0.0023148402339379995,
    #'Bax_dimerization_kf': 0.00087114734334575424,
    #'Bax_dimerization_kr': 0.0012875415792703037,
}

sc = 10

# Compartment-based vs. ODE
b1c = one_cpt.Builder(params_dict=params_dict)
b1c.build_model_tap1()
#b1c.build_model_ta()
fid = 1
b1c.run_model(tmax=tmax, figure_ids=[fid])

#bsc = site_cpt.Builder(scaling_factor=sc, params_dict=params_dict)
#bsc.build_model_ta()
#fid = 1
#bsc.run_model(tmax=tmax, figure_ids=[fid], num_sims=3)

bnc = n_cpt.Builder(scaling_factor=sc, params_dict=params_dict)
bnc.build_model_tap1()
#bnc.build_model_ta()
fid = 1
bnc.run_model(tmax=tmax, num_sims=12, figure_ids=[fid])
title('ODE vs. Compartment-based')

#legend()

# Site-based vs. ODE (Kappa)
#fid = 2 
#m1c.run_model(tmax=tmax, figure_ids=[fid])

#msc = tBid_Bax_sitec(scaling_factor=sc, params_dict=params_dict)
#msc.build_model0()
#mscx = msc.run_model(tmax=tmax, num_sims=2, use_kappa=True, figure_ids=[fid])
#title('ODE vs. Site-based (Kappa)')
#legend()

# Site-based vs. ODE (BNG)
#m1c.run_model(tmax=tmax, figure_ids=[1])

#msc = tBid_Bax_sitec(scaling_factor=10, params_dict=params_dict)
#msc.build_model0()
#msc.run_model(tmax=tmax, num_sims=3, use_kappa=False, figure_ids=[1])
#title('ODE vs. Site-based (BNG)')



