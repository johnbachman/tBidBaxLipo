from generate_model_ensemble_fit_files import generate_files

args = {'model': {
            'builder': ['one_cpt', 'lipo_sites'],
            'baxtranslocation': [1],
            'activation': [1],
            'reversal': [0, 1],
            'autoactivation': [0, 1],
            'dimerization': [0, 1, 2],
            'nbd': [1, 2],
            'bleach': [0],
            }
        }

def test_generate_files():
    generate_files(args)


