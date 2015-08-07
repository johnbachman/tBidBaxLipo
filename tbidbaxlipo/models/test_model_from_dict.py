from tbidbaxlipo.models.one_cpt import Builder
from nose.tools import *


@raises(KeyError)
def test_key_error_on_unknown_feature():
    model_dict = {'asdf': 1}
    bd = Builder()
    bd.build_model_from_dict(model_dict)

def test_model_name_sort_order():
    model_dict = {'activation': 1,
                  'baxtranslocation': 1}
    bd = Builder()
    bd.build_model_from_dict(model_dict)
    eq_(bd.model.name, 'Baxtr1Activ1')

def test_one_cpt_builder_abbrev():
    model_dict = {'builder': 'one_cpt',
                  'activation': 1,
                  'baxtranslocation': 1}
    bd = Builder()
    bd.build_model_from_dict(model_dict)
    eq_(bd.model.name, '1c_Baxtr1Activ1')

def test_lipo_sites_builder_abbrev():
    model_dict = {'builder': 'lipo_sites',
                  'activation': 1,
                  'baxtranslocation': 1}
    bd = Builder()
    bd.build_model_from_dict(model_dict)
    eq_(bd.model.name, 'ls_Baxtr1Activ1')

@raises(ValueError)
def test_error_on_unknown_implementation():
    model_dict = {'baxtranslocation': 9}
    bd = Builder()
    bd.build_model_from_dict(model_dict)

"""Here we check if the function is calling the right macros by checking
for the presence of named rules. Brittle if rule names later change."""
def test_translocation_1():
    """Check translocation by checking for Bax_mono_translocates_sol_to_ves."""
    model_dict = {'baxtranslocation': 1,
                  'bidtranslocation': 1}
    bd = Builder()
    bd.build_model_from_dict(model_dict)
    print bd.model.rules
    ok_(bd.model.rules['tBid_translocates_sol_to_ves'])
    ok_(bd.model.rules['Bax_mono_translocates_sol_to_ves'])

@raises(KeyError)
def test_error_on_bleach_without_nbd():
    model_dict = {'baxtranslocation': 1,
                  'activation': 1,
                  'bleach': 1}
    bd = Builder()
    bd.build_model_from_dict(model_dict)

