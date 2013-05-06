__author__ = "johnbachman"

from pylab import figure, title, legend
from tBid_Bax_1c import tBid_Bax_1c
from tBid_Bax_sitec import tBid_Bax_sitec
from tBid_Bax_nc import tBid_Bax_nc

tmax = 500

params_dict = {'Vesicles_0':1, 'tBid_0':1, 'Bax_0':50,
'Bax_transloc_kr':1e-2}

sc = 50

m1c = tBid_Bax_1c(params_dict=params_dict)
m1c.build_model0()

# Compartment-based vs. ODE
fid = 1
m1c.run_model(tmax=tmax, figure_ids=[fid])

#mnc = tBid_Bax_nc(scaling_factor=sc, params_dict=params_dict)
#mnc.build_model0()
#mncx = mnc.run_model(tmax=tmax, num_sims=3, figure_ids=[fid])
#title('ODE vs. Compartment-based')

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



