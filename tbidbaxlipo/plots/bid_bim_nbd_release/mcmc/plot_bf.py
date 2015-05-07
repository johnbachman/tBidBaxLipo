import glob
import sys
import pickle
import re
from itertools import product
from matplotlib import pyplot as plt
import numpy as np

# Store the bayes factor results as a dict of dicts
class EvidenceResults(object):
    def __init__(self):
        self.res_dict = {}

    def add(self, mutant, rep, num_confs, log_evidence):
        mut_dict = self.res_dict.setdefault(mutant, {})
        rep_dict = mut_dict.setdefault(rep, {})
        rep_dict[num_confs] = log_evidence

    def get(self, mutant, rep, num_confs):
        return self.res_dict[mutant][rep][num_confs]

    def get_result_arrays(self, key_list):
        evidence_list = []
        error_list = []
        for key_tuple in key_list:
            (evidence, error) = self.get(*key_tuple)
            evidence_list.append(evidence)
            error_list.append(error)
        return (np.array(evidence_list), np.array(error_list))

def load_files(files):
    er = EvidenceResults()
    for filename in files:
        print "Loading %s" % filename
        with open(filename) as f:
            (gf, sampler) = pickle.load(f)
        m = re.search(
                r'^pt_data1_Bid_NBD_([0-9]+)_r([1-3])_([2-5])confs_1.mcmc',
                filename)
        (mutant, rep, num_confs) = m.groups()
        log_ev = sampler.thermodynamic_integration_log_evidence()
        er.add(mutant, rep, num_confs, log_ev)
    return er

if __name__ == '__main__':
    #files = glob.glob('pt_data1_Bid_NBD_126_r1_2confs.mcmc')
    #files = [
    #    'pt_data1_Bid_NBD_126_r1_2confs_1.mcmc',
    #    'pt_data1_Bid_NBD_126_r1_3confs_1.mcmc',
    #]

    mut_list = ['54', '126']
    rep_list = ['2', '3']
    #rep_list = ['1']
    conf_list = ['3', '4', '5']
    #conf_list = ['2', '3', '4', '5']
    #conf_list = ['2', '3']
    key_list = [k for k in product(mut_list, rep_list, conf_list)]

    file_list = ['pt_data1_Bid_NBD_%s_r%s_%sconfs_1.mcmc' % subs
                 for subs in key_list]
    er = load_files(file_list)
    (ev, err) = er.get_result_arrays(key_list)
    print ev
    print err

    plt.ion()
    plt.figure()
    plt.bar(range(len(key_list)), -ev, yerr=err)
    plt.ylabel('-log(P(D|M))')

