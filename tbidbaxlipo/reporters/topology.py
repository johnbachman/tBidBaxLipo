import numpy as np
from bayessb.report import reporter, Result

reporter_group_name = 'Model topology'

def write_species(mcmc_set_name, model):
    """Write the species list of the model to a file.

    If the file does not exist, the function does nothing, but returns the name
    of the file.

    Parameters
    ----------
    mcmc_set_name : string
        The name of the MCMC set, which serves as the basename of the file to
        be written out.
    model : pysb.core.Model
        The model whose species are to be written to a file.

    Returns
    -------
    string
        The name of the species list filename.
    """
    species_filename = '%s_species.txt' % mcmc_set_name
    try:
        with open(species_filename, 'w') as f:
            #f.write('\n'.join([str(s) for s in model.rules]))
            #f.write('\n')
            f.write('\n'.join([str(s) for s in model.species]))
    except IOError as e:
        pass
    return species_filename

class TopologyResult(Result):
    """Implements specific HTML formatting for reporters indicating model
    topology."""
    def get_html(self):
        # Color-code the presence/absence of topology elements
        if self.value == True:
            color = "red"
        else:
            color = "lightgray"
        result_str = self.value

        # Add a link to the species list
        if self.link is not None:
            result_str = '<a href="%s">%s</a>' % (self.link, result_str)

        return '<td style="background: %s">%s</td>' % (color, result_str)

@reporter('Inserted Bax')
def inserted_Bax(mcmc_set):
    # Get model from first chain in the set
    model = mcmc_set.chains[0].options.model
    observable_names = [o.name for o in model.observables]

    # Write species list to a file
    species_filename = write_species(mcmc_set.name, model)

    if 'iBax' not in observable_names:
        return TopologyResult(False, species_filename)
    # This will be an empty list if the observable never occurs
    if model.observables['iBax'].species:
        return TopologyResult(True, species_filename)
    else:
        return TopologyResult(False, species_filename)

@reporter('Bax reverses')
def Bax_reverses(mcmc_set):
    # Get model from first chain in the set
    model = mcmc_set.chains[0].options.model

    # Write species list to a file
    species_filename = write_species(mcmc_set.name, model)

    if 'iBax_reverses' in [r.name for r in model.rules]:
        return TopologyResult(True, species_filename)
    else:
        return TopologyResult(False, species_filename)

@reporter('Bax/tBid inhibition ')
def tBid_inhibits_Bax(mcmc_set):
    # Get model from first chain in the set
    model = mcmc_set.chains[0].options.model

    # Write species list to a file
    species_filename = write_species(mcmc_set.name, model)

    if 'tBid_iBax_bind_at_bh3' in [r.name for r in model.rules]:
        return TopologyResult(True, species_filename)
    else:
        return TopologyResult(False, species_filename)

@reporter('Bax dimerizes')
def Bax_dimerizes(mcmc_set):
    # Get model from first chain in the set
    model = mcmc_set.chains[0].options.model
    observable_names = [o.name for o in model.observables]

    # Write species list to a file
    species_filename = write_species(mcmc_set.name, model)

    if 'Bax2' not in observable_names:
        return TopologyResult(False, species_filename)
    # This will be an empty list if the observable never occurs
    if model.observables['Bax2'].species:
        return TopologyResult(True, species_filename)
    else:
        return TopologyResult(False, species_filename)

@reporter('Bax tetramerizes')
def Bax_tetramerizes(mcmc_set):
    # Get model from first chain in the set
    model = mcmc_set.chains[0].options.model
    observable_names = [o.name for o in model.observables]

    # Write species list to a file
    species_filename = write_species(mcmc_set.name, model)

    if 'Bax4' not in observable_names:
        return TopologyResult(False, species_filename)
    # This will be an empty list if the observable never occurs
    if model.observables['Bax4'].species:
        return TopologyResult(True, species_filename)
    else:
        return TopologyResult(False, species_filename)

