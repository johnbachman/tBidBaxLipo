from nbd_mcmc_pysb import model_names #, nbd_site_names
from tbidbaxlipo.models.one_cpt import Builder
import subprocess

b = Builder()

nsteps = 300000

#model_names = ['tadt', 'tai']
nbd_site_names = ['c62']
# nbd_observables = [o.name for o in b.model.observables]
nbd_observables = ['Baxbh3', 'Bax2', 'tBidBax']

cmd_list = []
for model in model_names:
    for nbd_site_name in nbd_site_names:
        for nbd_observable in nbd_observables:
            cmd_list = ["python", "nbd_mcmc_pysb_lsf.py",
                        "model=%s" % model,
                        "nbd_site=%s" % nbd_site_name,
                        "nbd_observable=%s" % nbd_observable,
                        "nsteps=%d" % nsteps]
            print ' '.join(cmd_list)
            subprocess.call(cmd_list)

