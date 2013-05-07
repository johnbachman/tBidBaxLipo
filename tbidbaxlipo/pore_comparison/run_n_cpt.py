from pylab import figure, title, legend
from tbidbaxlipo.models import one_cpt, site_cpt, n_cpt
from pysb import bng

tmax = 500
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

bnc = n_cpt.Builder(scaling_factor=sc, params_dict=params_dict)
bnc.build_model_tap1()
bng.run_ssa(bnc.model, t_end=tmax, n_steps=100, cleanup=False, output_dir='.')

#fid = 1
#bnc.run_model(tmax=tmax, num_sims=12, figure_ids=[fid])



