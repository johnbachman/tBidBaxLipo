#from nbd_mcmc_pysb import model_names
from tbidbaxlipo.models.one_cpt import Builder
import subprocess

model_names = ['tairdt']

b = Builder()

nsteps = 150000

nbd_site_combos = [('c3', 'c62')]
nbd_obs_combos = [('iBax', 'Baxbh3'), ('iBax', 'Bax2')]

cmd_list = []
for model in model_names:
    for nbd_site_combo in nbd_site_combos:
        nbd_site_str = '-'.join(nbd_site_combo)
        for nbd_obs_combo in nbd_obs_combos:
            nbd_obs_str = '-'.join(nbd_obs_combo)
            cmd_list = ["python", "-m", "tbidbaxlipo.nbd_mcmc_pysb_lsf",
                        "model=%s" % model,
                        "nbd_sites=%s" % nbd_site_str,
                        "nbd_observables=%s" % nbd_obs_str,
                        "nsteps=%d" % nsteps]
            print ' '.join(cmd_list)
            subprocess.call(cmd_list)

