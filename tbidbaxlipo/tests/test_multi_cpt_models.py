from tbidbaxlipo.models.multi_cpt import Builder
from nose.tools import ok_

def test_build_model_t():
    b = Builder()
    b.build_model_t()
    b.make_multi_compartment()
    ok_(check_monomer_and_parameter_refs(b.model),
        "Monomer reference check failed.")

def check_monomer_and_parameter_refs(model):
    # Check rules
    for rule in model.rules:
        # For both the left and right hand sides, get the complex patterns...
        lh_cps = rule.rule_expression.reactant_pattern.complex_patterns
        rh_cps = rule.rule_expression.product_pattern.complex_patterns
        for side_cps in (lh_cps, rh_cps):
            # For all the complex patterns on this side of the rule...
            for cp in side_cps:
                # For all the monomer patterns in this complex pattern...
                for mp in cp.monomer_patterns:
                    # Check that the monomer in the monomer pattern is in
                    # the model's list of monomers
                    if mp.monomer not in model.monomers:
                        print("Monomer %s in rule %s not in "
                              "model.monomers" % (mp.monomer, rule))
                        return False
        # Check that the parameter(s) are in model.parameters
        if rule.rate_forward is not None and \
           rule.rate_forward not in model.parameters:
            print("Forward rate parameter %s for rule %s not in "
                  "model.parameters" % (rule.rate_forward, rule))
            return False
        if rule.rate_reverse is not None and \
           rule.rate_reverse not in model.parameters:
            print("Reverse rate parameter %s for rule %s not in "
                  "model.parameters" % (rule.rate_forward, rule))
            return False

    # Check observables
    for obs in model.observables:
        ob_cps = obs.reaction_pattern.complex_patterns
        for cp in ob_cps:
            for mp in cp.monomer_patterns:
                if mp.monomer not in model.monomers:
                    print("Monomer %s in observable %s not in model.monomers"
                          % (mp.monomer, obs))
                    return False

    # Check initial conditions
    for (ic_cp, ic_value) in model.initial_conditions:
        for mp in ic_cp.monomer_patterns:
            if mp.monomer not in model.monomers:
                print("Monomer %s in initial condition (%s, %s) not in "
                      "model.monomers" % (mp.monomer, ic_cp, ic_value))
                return False

    # We made it this far, so all of our references have checked out!
    return True
