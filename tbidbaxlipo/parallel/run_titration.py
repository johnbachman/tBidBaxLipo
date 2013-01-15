import numpy as np

conditions = [{'tBid_0': 10, 'Bax_0': 30, 'Vesicles_0': 50},
              {'tBid_0': 20, 'Bax_0': 100, 'Vesicles_0': 5}]

num_sims = 4

def submit_jobs():
    for (cond_num, conditions_dict) in enumerate(conditions):
        for sim_num in range(0, num_sims):
            print('bsub python run_site_cpt.py -basename %s_cond%d_sim%d'
                   % ('basename', cond_num, sim_num))

