class ParameterSet(object):
    # should be able to load a parameter set for a model easily
    # model 3r
    # a list of parameters, or a dict with
    # the parameter names and their values
    #[  1.16737842e+02   6.54218946e+02   7.22384624e+05 
    # need a way to specify which parameters should be set and which
    # not, presumably done by the obj_func; if the icparams
    # are not to be fit, then obj_func should not return them
    # along in the param dict
    def __init__(self, param_dict, best_err, model=None):
        self.param_dict = param_dict
        self.best_err = best_err
        self.model = model

"""
Jmin:  0.0863839628828
T:  1751
feval:  34
iters:  369
"""
# Would also be nice if the function could interpret the changes
# in the parameters resulting from fitting, i.e., a better fit
# was achieved by making binding weaker or on-rates slower, etc.

ps1 = ParameterSet({'Bax_transloc_kf': 3.45078952e-04,
             'Bax_transloc_kr': 1.16898702e+00,
             'tBid_transloc_kf': 1.81447359e-02,
             'tBid_transloc_kr': 2.13787673e+00,
             'tBid_mBax_kf': 6.82463546e-03,
             'tBid_mBax_kr': 2.69976698e-01,
             'mBaxtBid_to_iBaxtBid_k': 1.44906963e-03,
             'tBid_iBax_kr': 6.54536144e-02,
             'tBid_iBax_kf': 1.11620976e-02,
             'iBaxtBid_to_mBaxtBid_k': 5.49437163e-03,
             'Bax_dimerization_kf': 1.38217318e-02,
             'k_eflx': 1.08106987e+04,
             'lipo_eflx': 1.73086744e+02,
             'Bax_dimerization_kr': 6.30226287e-01},
              58.1016586865)
# min obj func val = 
#init_conds = {'tBid_0': 15, 'Vesicles_0': 10}


