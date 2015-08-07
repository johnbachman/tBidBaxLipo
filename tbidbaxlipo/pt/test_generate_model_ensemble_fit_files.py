from generate_model_ensemble_fit_files import generate_files
from copy import copy

args_no_builder = {
    'model': {
        'baxtranslocation': [1],
        'activation': [1],
        'reversal': [0, 1],
        'autoactivation': [0, 1],
        'dimerization': [0, 1, 2],
        'nbd': [1, 2],
        'bleach': [0]}}

args_with_builder = copy(args_no_builder)
args_with_builder['model']['builder'] = ['one_cpt', 'lipo_sites']

def test_generate_files_no_builder():
    generate_files(args_no_builder, 'test')

def test_generate_files_with_builder():
    generate_files(args_with_builder, 'test')


