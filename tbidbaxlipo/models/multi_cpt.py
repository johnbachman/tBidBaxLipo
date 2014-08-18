from tbidbaxlipo.models import core, one_cpt
from bayessb.priors import Normal
from pysb.core import Model, Monomer, ComponentSet
from copy import copy, deepcopy

class Builder(core.Builder):
    def make_multi_compartment(self):
        # Multi-compartment model
        mc_model = Model(base=self.model)

        # MONOMERS
        # First we need to iterate over the monomers and figure out which
        # ones can be localized to different vesicles (determined by having
        # a site 'cpt' with a possible value 'ves'). If the Monomer has this
        # property, we copy it to the new model but specify that it can take
        # on n values, for the n possible vesicle compartments.
        mc_monomers = ComponentSet()
        monomer_map = {}
        for m in mc_model.monomers:
            if 'cpt' in m.sites and \
               'ves' in m.site_states['cpt']:
                mc_monomer = deepcopy(m)
                # FIXME Need to make iterative, to cover all of the vesicles
                mc_monomer.site_states['cpt'] = ['ves1']
                mc_monomers.add(mc_monomer)
                monomer_map[m] = mc_monomer
            else:
                mc_monomers.add(m)
                monomer_map[m] = m
        # Overwrite the Monomer ComponentSet in the copied model with the newly
        # updated set
        mc_model.monomers = mc_monomers

        # RULES
        mc_rules = ComponentSet()
        for rule in mc_model.rules:
            # Need to do a copy, rather than a deep copy, because otherwise we
            # replace all the monomer references with new copies that won't
            # sync with our monomer list
            mc_rule = copy(rule)
            lh_cps = mc_rule.rule_expression.reactant_pattern.complex_patterns
            rh_cps = mc_rule.rule_expression.product_pattern.complex_patterns
            # For both the left and right hand sides, get the complex patterns
            for side_cps in (lh_cps, rh_cps):
                # For all the complex patterns on this side of the rule...
                for cp in side_cps:
                    # For all the monomer patterns in this complex pattern...
                    for mp in cp.monomer_patterns:
                        # Replace references to 'cpt=ves' with multiple
                        # compartment references
                        if 'cpt' in mp.site_conditions and \
                           mp.site_conditions['cpt'] == 'ves':
                            # FIXME Need to make this iterative, to cover all
                            # of the vesicles
                            mp.site_conditions['cpt'] = 'ves1'
                        # Substitute the appropriate monomer reference, as
                        # given by the monomer_map, initialized above
                        mp.monomer = monomer_map[mp.monomer]

            # Add the rule to the new rule component set
            mc_rules.add(mc_rule)
        # Reinitialize the rule component set of the model with our updated one
        mc_model.rules = mc_rules

        # INITIAL CONDITIONS
        mc_ic_list = []
        for (ic_cp, ic_value) in mc_model.initial_conditions:
            mc_ic_cp = copy(ic_cp)
            for mp in mc_ic_cp.monomer_patterns:
                # Replace references to 'cpt=ves' with multiple
                # compartment references
                # FIXME Need to make this iterative, to cover all
                # of the vesicles
                if 'cpt' in mp.site_conditions and \
                   mp.site_conditions['cpt'] == 'ves':
                    mp.site_conditions['cpt'] = 'ves1'
                # Substitute the appropriate monomer reference, as
                # given by the monomer_map, initialized above
                mp.monomer = monomer_map[mp.monomer]
            mc_ic_list.append((mc_ic_cp, ic_value))
        # Reinitialize the initial condition list with the updated one
        mc_model.initial_conditions = mc_ic_list

        # OBSERVABLES
        mc_obs_set = ComponentSet()
        for obs in mc_model.observables:
            # Again, don't do deep copy (see note above for rules)
            mc_obs = copy(obs)
            ob_cps = mc_obs.reaction_pattern.complex_patterns
            for cp in ob_cps:
                for mp in cp.monomer_patterns:
                    if 'cpt' in mp.site_conditions and \
                       mp.site_conditions['cpt'] == 'ves':
                        mp.site_conditions['cpt'] = 'ves1'
                    # Substitute the appropriate monomer reference, as
                    # given by the monomer_map, initialized above
                    mp.monomer = monomer_map[mp.monomer]
            # Add the updated observable to the new observable component set
            mc_obs_set.add(mc_obs)
        # Reinitialize the observables of the model to the new set
        mc_model.observables = mc_obs_set

        # Reset the builder's model object to the multi-compartmentalized model
        self.model = mc_model

if __name__ == '__main__':
    from pysb.integrate import Solver
    from matplotlib import pyplot as plt
    import numpy as np

    b = Builder()
    b.build_model_t()
    b.make_multi_compartment()
    """
    t = np.linspace(0, 1e5)
    s = Solver(b.model, t)
    s.run()
    plt.ion()
    plt.figure()
    plt.plot(t, s.yexpr['NBD'])
    """
