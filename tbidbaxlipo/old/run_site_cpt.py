"""Code to run a single simulation of the desired model at the command line.
Since it can be called at the command line it can be parallelized easily
by submitting multiple simultaneous jobs through bsub.
"""

__author__ = "johnbachman"

from tbidbaxlipo.models import site_cpt


#params_dict = {'Vesicles_0':1, 'tBid_0':1, 'Bax_0':50,
#'Bax_transloc_kr':1e-2}

# Make the model we want to run
def run(output_dir, base_filename):
    tmax = 1000
    mb = site_cpt.Builder()
    mb.build_model0()
    mb.run_model(tmax=tmax, base_filename=base_filename)

if __name__ == '__main__':
    run()


