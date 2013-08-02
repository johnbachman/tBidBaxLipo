from pore_plots import *
from tbidbaxlipo.models import one_cpt

# No tBid/iBax binding
params_dict = {'tBid_mBax_kf':1e-1, 'tBid_iBax_kc':1e-2}
b = one_cpt.Builder(params_dict=params_dict)

b.translocate_tBid_Bax()
b.tBid_activates_Bax()
b.Bax_dimerizes()
b.pores_from_Bax_dimers()

plot_bax_titration(b.model)

b.model.parameters['Bax_0'].value = 100
plot_tbid_titration(b.model)

# With tBid/iBax binding
b = one_cpt.Builder(params_dict=params_dict)

b.translocate_tBid_Bax()
b.tBid_activates_Bax()
b.iBax_binds_tBid_at_bh3()
b.Bax_dimerizes()
b.pores_from_Bax_dimers()

plot_bax_titration(b.model)

b.model.parameters['Bax_0'].value = 100
plot_tbid_titration(b.model)

