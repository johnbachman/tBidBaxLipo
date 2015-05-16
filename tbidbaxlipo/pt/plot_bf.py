import glob
import sys
import pickle
import re
from itertools import product
from matplotlib import pyplot as plt
import numpy as np
import csv

# Store the bayes factor results as a dict of dicts
class EvidenceResults(object):
    """Stores the table of evidence results from parallel tempering runs.

    Attributes
    ----------
    res_dict : dict of dicts of dicts
        The top level dict is keyed by activator, then mutant, then replicate,
        then the number of conformations in the model. At the bottom level
        the dict contains tuples of (log evidence, error) as returned by
        emcee.PTSampler.thermodynamic_integration_log_evidence().
    """

    def __init__(self):
        self.res_dict = {}

    def add(self, activator, mutant, rep, num_confs, log_evidence):
        """Add an evidence result for the activator/mutant/rep/model."""
        act_dict = self.res_dict.setdefault(activator, {})
        mut_dict = act_dict.setdefault(mutant, {})
        rep_dict = mut_dict.setdefault(rep, {})
        rep_dict[num_confs] = log_evidence

    def get(self, activator, mutant, rep, num_confs):
        """Get the evidence result for the activator/mutant/rep/model."""
        return self.res_dict[activator][mutant][rep][num_confs]

    def get_result_arrays(self, key_list):
        """Get a list of evidence results for a given set of keys.

        Parameters
        ----------
        key_list : list of tuples
            Each tuple in the list should contain an activator, mutant, rep,
            and number of conformations, to be passed as arguments to
            :py:meth:`EvidenceResults.get`.

        Returns
        -------
        tuple : (numpy.array, numpy.array)
            The first numpy.array in the tuple contains the log evidence
            values; the second contains the corresponding errors.
        """

        evidence_list = []
        error_list = []
        for key_tuple in key_list:
            (evidence, error) = self.get(*key_tuple)
            evidence_list.append(evidence)
            error_list.append(error)
        return (np.array(evidence_list), np.array(error_list))

    def get_results_table(self, conf_list):
        """Generate a table of all results as a list of lists."""
        # Iterate over the activators (sorted)
        rows = []
        for act in sorted(self.res_dict.keys()):
            act_dict = self.res_dict[act]
            for mut in sorted(act_dict.keys(), key=lambda mut: int(mut)):
                mut_dict = act_dict[mut]
                for rep in sorted(mut_dict.keys(), key=lambda rep: int(rep)):
                    rep_dict = mut_dict[rep]
                    row = [act, mut, rep]
                    for conf in conf_list:
                        try:
                            (log_ev, err) = rep_dict[conf]
                            row.append(log_ev)
                            row.append(err)
                        # If we don't yet have results for this conformation,
                        # put None in the table
                        except KeyError:
                            row.append(None)
                            row.append(None)
                    rows.append(row)
        return rows

    def write_csv(self, filename, conf_list):
        results = self.get_results_table(conf_list)
        with open(filename, 'w') as csvfile:
            cw = csv.writer(csvfile, delimiter=',')
            for row in results:
                cw.writerow(row)
        return

def load_files(files):
    """For a given list of .mcmc files, load and store the evidence results."""
    er = EvidenceResults()
    for filename in files:
        print "Loading %s" % filename
        try:
            with open(filename) as f:
                (gf, sampler) = pickle.load(f)
            m = re.search(
                 r'^pt_data1_([Bidm]+)_NBD_([0-9]+)_r([1-3])_([2-5])confs.mcmc',
                 filename)
            (activator, mutant, rep, num_confs) = m.groups()
            log_ev = sampler.thermodynamic_integration_log_evidence()
            er.add(activator, mutant, rep, num_confs, log_ev)
        except IOError:
            print("File %s not found, skipping" % filename)
    return er

if __name__ == '__main__':
    #files = glob.glob('pt_data1_Bid_NBD_126_r1_2confs.mcmc')
    #files = [
    #    'pt_data1_Bid_NBD_126_r1_2confs_1.mcmc',
    #    'pt_data1_Bid_NBD_126_r1_3confs_1.mcmc',
    #]
    from tbidbaxlipo.data.parse_bid_bim_nbd_release import nbd_residues
    # Filter out the WT residue
    mut_list = [res for res in nbd_residues if res != 'WT']
    rep_list = [str(n) for n in range(1, 4)]
    conf_list = [str(n) for n in range(2, 6)]
    key_list = [k for k in product(mut_list, rep_list, conf_list)]
    file_list = ['pt_data1_Bid_NBD_%s_r%s_%sconfs.mcmc' % subs
                 for subs in key_list]
    ev_results = load_files(file_list)
    ev_results.write_csv('results.csv', conf_list)
    #(ev, err) = er.get_result_arrays(key_list)
    #print ev
    #print err

    #plt.ion()
    #plt.figure()
    #plt.bar(range(len(key_list)), -ev, yerr=err)
    #plt.ylabel('-log(P(D|M))')

