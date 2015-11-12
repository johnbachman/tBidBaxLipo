import cPickle
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import \
            set_fig_params_for_publication, format_axis, fontsize
import re
from tbidbaxlipo.models import model_aliases

set_fig_params_for_publication()

with open('evidence_list.pck') as f:
    ev_list = cPickle.load(f)

# Group 1:
group1 = [
    'pt_140318_1c_Baxtr1Activ1Rever1Nbd1.mcmc',
    'pt_140318_1c_Baxtr1Activ1Nbd1.mcmc',
    'pt_140318_ls_Baxtr1Activ1Rever1Nbd1.mcmc',
    'pt_140318_ls_Baxtr1Activ1Nbd1.mcmc']

group2 = [
    'pt_140318_1c_Baxtr1Activ1Nbd1.mcmc',
    'pt_140318_1c_Baxtr1Activ1Rever1Nbd1.mcmc',
    'pt_140318_1c_Baxtr1Activ1Dimer2Nbd2.mcmc',
    'pt_140318_1c_Baxtr1Activ1Dimer1Nbd3.mcmc',
    'pt_140318_1c_Baxtr1Activ1Rever1Dimer2Nbd2.mcmc',
    'pt_140318_1c_Baxtr1Activ1Dimer2Nbd3.mcmc',
    'pt_140318_1c_Baxtr1Activ1Rever1Dimer1Nbd2.mcmc',
    'pt_140318_1c_Baxtr1Activ1Rever1Dimer1Nbd3.mcmc',
    'pt_140318_1c_Baxtr1Activ1Dimer1Nbd2.mcmc',
    'pt_140318_1c_Baxtr1Activ1Rever1Dimer2Nbd3.mcmc',
    'pt_140318_1c_Baxtr1Activ1Rever1Dimer1Nbd1.mcmc',
    'pt_140318_1c_Baxtr1Activ1Rever1Dimer2Nbd1.mcmc',
    'pt_140318_1c_Baxtr1Activ1Dimer2Nbd1.mcmc',
    'pt_140318_1c_Baxtr1Activ1Dimer1Nbd1.mcmc',]

ev_list1 = [ev_list_item for ev_list_item in ev_list
               if ev_list_item[0] in group1]

names, ev = zip(*ev_list1)
ev_vals, ev_errs = zip(*ev)
ev_vals = np.max(ev_vals) - np.array(ev_vals)

plt.ion()
plt.figure(figsize=(1.5, 1), dpi=300)

"""
# Group 1: Dimerization, NBD species is dimer
g1_ix = 0
g2_ix = 8
g3_ix = 12

# Group 1: Dimerization, NBD species is dimer
plt.barh(range(g2_ix), ev_vals[0:g2_ix], color='r')
# Group 2: No dimerization or NBD species is monomer
plt.barh(range(g2_ix, g3_ix), ev_vals[g2_ix:g3_ix], color='g')
# Group 3: Liposome binding sites
plt.barh(12, ev_vals[12], color='b')
plt.barh(14, ev_vals[14], color='b')
# Group 4: ???
plt.barh(13, ev_vals[13], color='g')
plt.barh(15, ev_vals[15], color='g')
"""
plt.barh(range(len(ev_vals)), ev_vals, xerr=ev_errs, ecolor='k')
plt.xlabel(r'$\Delta$ log(Marg. Lkl.)')
plt.xticks(rotation='vertical')
parsed_names = [model_aliases[re.match('pt_140318_(\w*)\.mcmc', name).groups()[0]]
                for name in names]
plt.yticks(np.arange(len(ev_vals)) + 0.4, parsed_names)
plt.xlim(left=0)
plt.subplots_adjust(bottom=0.42, left=0.17, right=0.95, top=0.95)
ax = plt.gca()
format_axis(ax)

plt.savefig('140318_evidence_barplot1.pdf', dpi=300)
plt.savefig('140318_evidence_barplot1.png', dpi=300)

# Make figure 2
ev_list2 = [ev_list_item for ev_list_item in ev_list
               if ev_list_item[0] in group2]

names, ev = zip(*ev_list2)
ev_vals, ev_errs = zip(*ev)
ev_vals = np.max(ev_vals) - np.array(ev_vals)

plt.figure(figsize=(1.5, 3), dpi=300)

"""
# Group 1: Dimerization, NBD species is dimer
g1_ix = 0
g2_ix = 8
g3_ix = 12

# Group 1: Dimerization, NBD species is dimer
plt.barh(range(g2_ix), ev_vals[0:g2_ix], color='r')
# Group 2: No dimerization or NBD species is monomer
plt.barh(range(g2_ix, g3_ix), ev_vals[g2_ix:g3_ix], color='g')
# Group 3: Liposome binding sites
plt.barh(12, ev_vals[12], color='b')
plt.barh(14, ev_vals[14], color='b')
# Group 4: ???
plt.barh(13, ev_vals[13], color='g')
plt.barh(15, ev_vals[15], color='g')
"""
plt.barh(range(len(ev_vals)), ev_vals, xerr=ev_errs, ecolor='k')
plt.xlabel(r'$\Delta$ log(Marg. Lkl.)')
plt.xticks(rotation='vertical')
parsed_names = [model_aliases[re.match('pt_140318_(\w*)\.mcmc', name).groups()[0]]
                for name in names]
plt.xlim(left=0)
plt.yticks(np.arange(len(ev_vals)) + 0.4, parsed_names)
plt.subplots_adjust(bottom=0.21, left=0.17, right=0.95, top=0.95)
ax = plt.gca()
format_axis(ax)

plt.savefig('140318_evidence_barplot2.pdf', dpi=300)
plt.savefig('140318_evidence_barplot2.png', dpi=300)

