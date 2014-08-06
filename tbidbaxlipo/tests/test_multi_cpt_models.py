from tbidbaxlipo.models.multi_cpt import Builder
from nose.tools import ok_

def test_build_model_t():
    b = Builder()
    b.build_model_t()
    b.make_multi_compartment()
    ok_(check_monomer_references())

def check_monomer_references():
    return False
