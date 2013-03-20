from nbd_mcmc_pysb import model_names
from tbidbaxlipo.models.one_cpt import Builder
import subprocess

b = Builder()

nsteps = 150000

nbd_site_names = ['c62']
nbd_observables = ['iBax', 'Baxbh3', 'Bax2', 'tBidBax']

cmd_list = []
for model in model_names:
    for nbd_site_name in nbd_site_names:
        for nbd_observable in nbd_observables:
            cmd_list = ["python", "-m", "tbidbaxlipo.nbd_mcmc_pysb_lsf",
                        "model=%s" % model,
                        "nbd_site=%s" % nbd_site_name,
                        "nbd_observable=%s" % nbd_observable,
                        "nsteps=%d" % nsteps]
            print ' '.join(cmd_list)
            subprocess.call(cmd_list)

