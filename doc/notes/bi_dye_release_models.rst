Bimolecular (tBid/Bax) dye release models
=========================================

one_cpt and n_cpt agreement:

* Interestingly, with pores from Bax monomers, the one_cpt and n_cpt models
  agree fairly closely when there are more tBid molecules than vesicles, EVEN
  when the off-rate is 0.
* They don't agree when the concentration of tBid is equal to that of vesicles,
  AND the off-rate is zero. In these cases the number of pores per vesicle
  is actually equal, but the dye release calculation diverges substantially.
  This is because a substantial number of pores are formed on vesicles that
  already have pores (violation of the Poisson assumption).

