from tbidbaxlipo.models.one_cpt import Builder
from nose.tools import *


@raises(KeyError)
def test_key_error_on_unknown_feature():
    model_dict = {'asdf': 1}
    bd = Builder()
    bd.build_model_from_dict(model_dict)

def test_model_name_sort_order():
    model_dict = {'activation': 1,
                  'translocation': 1}
    bd = Builder()
    bd.build_model_from_dict(model_dict)
    eq_(bd.model_name, 'Tra1Act1')

@raises(ValueError)
def test_error_on_unknown_implementation():
    model_dict = {'translocation': 9}
    bd = Builder()
    bd.build_model_from_dict(model_dict)

"""Here we check if the function is calling the right macros by checking
for the presence of named rules. Brittle if rule names later change."""
def test_translocation_1():
    """Check translocation by checking for Bax_mono_translocates_sol_to_ves."""
    model_dict = {'translocation': 1}
    bd = Builder()
    bd.build_model_from_dict(model_dict)
    ok_(bd.model.rules['tBid_translocates_sol_to_ves'])
    ok_(bd.model.rules['Bax_mono_translocates_sol_to_ves'])

@raises(KeyError)
def test_error_on_bleach_without_nbd():
    model_dict = {'translocation': 1,
                  'activation': 1,
                  'bleach': 1}
    bd = Builder()
    bd.build_model_from_dict(model_dict)

