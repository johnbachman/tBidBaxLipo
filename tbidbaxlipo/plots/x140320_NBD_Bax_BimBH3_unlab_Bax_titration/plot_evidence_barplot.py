import cPickle
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import \
            set_fig_params_for_publication, format_axis, fontsize
import re
from tbidbaxlipo.models import model_aliases
import sys

set_fig_params_for_publication()

evidence_file = sys.argv[1]
with open(evidence_file) as f:
    ev_list = cPickle.load(f)

names, ev = zip(*ev_list)
ev_vals, ev_errs = zip(*ev)
ev_vals = np.max(ev_vals) - np.array(ev_vals)

plt.ion()
plt.figure(figsize=(1.5, 3), dpi=300)

# Group 1: Dimerization, NBD species is dimer
#g1_ix = 0
#g2_ix = 8
#g3_ix = 12

# Group 1: Dimerization, NBD species is dimer
#plt.bar(range(g2_ix), ev_vals[0:g2_ix], color='r')
# Group 2: No dimerization or NBD species is monomer
#plt.bar(range(g2_ix, g3_ix), ev_vals[g2_ix:g3_ix], color='g')
# Group 3: Liposome binding sites
#plt.bar(12, ev_vals[12], color='b')
#plt.bar(14, ev_vals[14], color='b')
# Group 4: ???
#plt.bar(13, ev_vals[13], color='g')
#plt.bar(15, ev_vals[15], color='g')
parsed_names = [model_aliases[re.match('pt_140320_(\w*)\.mcmc', name).groups()[0]]
                for name in names]
for m_ix, name in enumerate(parsed_names):
    if name in ['M%d' % i for i in range(9, 13)]:
        plt.barh(m_ix, ev_vals[m_ix], color='r', xerr=ev_errs[m_ix],
                 ecolor='k')
    elif name in ['M%d' % i for i in range(13, 17)]:
        plt.barh(m_ix, ev_vals[m_ix], color='m', xerr=ev_errs[m_ix],
                 ecolor='k')
    elif name in ['M%d' % i for i in range(4, 9)]:
        plt.barh(m_ix, ev_vals[m_ix], color='g', xerr=ev_errs[m_ix],
                 ecolor='k')
    else:
        plt.barh(m_ix, ev_vals[m_ix], color='b', xerr=ev_errs[m_ix],
                 ecolor='k')

plt.xlabel(r'$\mathsf{\Delta}$ log(Marg. Lkl.)')
plt.xticks(rotation='vertical')
plt.xlim(left=0)
plt.yticks(np.arange(len(ev_vals)) + 0.4, parsed_names)
plt.subplots_adjust(bottom=0.21, left=0.17, right=0.95, top=0.95)
ax = plt.gca()
format_axis(ax)

plt.savefig('140320_evidence_barplot.pdf', dpi=300)
plt.savefig('140320_evidence_barplot.png', dpi=300)
