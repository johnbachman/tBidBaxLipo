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
    def __init__(self, param_dict, name=None, best_err=None, model=None):
        self.name = name
        self.param_dict = param_dict
        self.best_err = best_err
        self.model = model

    def diff(self, other):
        """ Compares the parameter values in this parameter set with those
            in the other parameter set. """
        # Iterate over the parameters in this set
        for i, param in enumerate(self.param_dict):
            if (param in other.param_dict):
                my_val = self.param_dict[param]
                other_val = other.param_dict[param]            
                output = "%s: %f --> %f" % (param, my_val, other_val)
                if (other_val < my_val):
                    output += "\033[0;31m slower\033[0m"
                elif (other_val > my_val):
                    output += "\033[0;32m faster\033[0m"
                print output
            else:
                print("WARNING: The parameter " + param + " is not contained in " +
                      "the other parameter set.")

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

# Parameter set resulting from fitting this model
"""def build_model1():
    print "Building model 1: Activation, dimerization"
    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax(bax_site='a6')
    Bax_dimerizes(dimer_diss_rate=2.5e-1)
    dye_release_dimeric()
"""
# to the dose response at 2, 21.5, and 40nM tBid.
ps_m1 = ParameterSet({
 'Bax_transloc_kf': 1.54490368e-03,
 'Bax_transloc_kr': 4.50607756e-02,
 'tBid_transloc_kf': 8.94212792e-03,
 'tBid_transloc_kr': 9.93132018e+00,
 'tBid_mBax_kf': 1.08401291e-04,
 'tBid_mBax_kr': 2.49910419e-02,
 'mBaxtBid_to_iBaxtBid_k': 8.50105160e-03,
 'tBid_iBax_kr': 3.85287413e-03,
 'tBid_iBax_kf': 1.65033482e-02,
 'Bax_dimerization_kf': 1.51859505e+00,
 'Bax_dimerization_kr': 3.05646033e-03,
 'k_eflx': 4.13273113e+03,
 'lipo_eflx': 6.15208340e+01},
 best_err=265.011127508)

# min obj func val = 
#init_conds = {'tBid_0': 15, 'Vesicles_0': 10}

"""What would you need to recreate an analysis?
- The model (e.g., m.model, with initial values)
- The fitting function (e.g., o.fit_grid())
- The data (e.g., g.data_bgsub)
- Any random number generator seeds
- Then, to avoid having to recreate, also want
   - The final parameter set
   - The val of the objective function
   - Any plots (before, after)
- Comments
"""
