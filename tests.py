"""
PySB tests meant to represent the knowledge about the dynamics of Bcl-2 family
interactions as studied using liposomal model systems.
"""
import pygraphviz as pgv
import networkx as nx

from pysb import kappa
from pysb.annotation import Citation, Context
from nose.tools import assert_true

from tBid_Bax_sitec import tBid_Bax_sitec

m = tBid_Bax_sitec(scaling_factor=1)
m.build_model0()

# Run the static analysis and get the graph back through PyGraphViz

# Set model conditions to those of the experiment
m['Vesicles_0'].value = 1
m['Bax_0'].value = 1

im_filename = kappa.influence_map(m.model)
inf_map = nx.Graph(pgv.AGraph(im_filename))

# Yethon et al. paper (2003)
# ==========================

yethon = Citation(
    """Interaction with a Membrane Surface Triggers a Reversible Conformational
    Change in Bax Normally Associated with Induction of Apoptosis""",
    """Jeremy A. Yethon, Raquel F. Epand, Brian Leber, Richard M. Epand,
    and David W. Andrews""",
    "14522999"
    )

# Observations from this paper, roughly in the order in which they are
# mentioned:

def test_lipids_trigger_6A7_exposure():
    """Transient association with lipids leads to 6A7 expos (3C)
    exposure, but not insertion, oligomerization, or permeabilization.
    """
    # Get a model
    assert_true(True)
    # Perform static test
    # Bax(loc='c', c6A7='n', bax2=None)
    #b = (nx.has_path(inf_map, Bax(c6A7='y')) and
    #    not nx.has_path(Bax(ins='y')) and
    #    not nx.has_path(Bax(bax2=1) % Bax(bax2=1)))

    # Perform dynamic test

# * Bax at the ER is involved in regulation of calcium ion fluxes (ref)
#
# * Ku70 interaction with Bax in the cytoplasm prevents its
#   translocation to membranes (ref)
#
# * Membrane permeabilization by Bax depends on intrinsic lipid curvature (ref)
#
# * Membrane permeabilization by Bax depends on cardiolipin (ref)
#
# * tBid-Bax pores can release dextrans up to 2MDa (ref)
#
# * The lipid-induced change is not induced solely by electrostatic interactions.
#   (not shown)
#
# * Bax activation can be inhibited by the membrane-active drug propanolol? (ref)
#   ...and is dependent on Mg2+ and Ca2+ in some systems? However, these were
#   not shown to affect Bax 6A7 (could test propanolol dependence in release
#   assay)
#
# * Neither mitochondria or ER membranes alone triggered conformational
#   change that was associated with liposomes. (Fig 1)
#
# * Negligible binding of Bax to liposomes in absence of tBid (fig 4?)

def no_binding_of_Bax_to_lipos_without_tBid():
    pass

# * In the absence of tBid, Bax doesn't insert or oligomerize (Fig 4).

def no_insertion_of_Bax_without_tBid():
    pass

def no_oligomerization_of_Bax_without_tBid():
    pass

# * A very small fraction of Bax incubated with liposomes forms dimers that
#   survive lysis of liposomes with CHAPS (interesting!). Are these dimers
#   inserted or not? (Fig 5A)

# * Bax alone with liposomes triggers minimal (<10%) ANTS release (Fig 5B).
#   Conditions: 2h incubation
def Bax_alone_triggers_minimal_dye_release_at_2hr():
    pass

# * tBid was only able to induce activation of Bax in the presence of liposomes
#   (Fig 6D)
def tBid_only_activates_Bax_with_lipos():
    pass

# * BclXL + liposomes does not block the 6A7 exposure (Fig 6B, lane 3).
#   One thing this doesn't rule out though, is that activated BclXL
#   (pre-equilibrated with tBid + liposomes, for example, doesn't prevent
#   6A7 exposure. Perhaps this is dealt with in the Billen paper.
def BclXL_does_not_prevent_6A7_exposure():
    pass

# * BclXL (200nM) + 20nM tBid + 100nM Bax almost completely prevents ANTS exposure
def BclXL_prevents_dye_release():
    pass

# = Discussion =
# * Contact with mitochondrial or ER membranes doesn't trigger Bax 6A7 exposure--
#   perhaps because of excessive protein blocking contact with the membranes?
#
# * Bax 6A7 reverses in cells after pro-apoptotic stimulus is removed.
#
# * "They conclude that the C-terminal Bax helix is a signal-anchor required
#   for membrane targeting and insertion but that elements at the N terminus are
#   key to regulating this process. Deletion of the N-terminal 20 amino acids of
#   Bax (encompassing the 6A7 epitope) stimulated Bax targeting to membranes in
#   vitro and increased the pro-apoptotic activity of Bax in vivo (4).
#
# * Found that 6A7 exposure didn't require particular lipid composition,
#   seemingly contradictory given results showing cardiolipin necessity.
#   "However, cardiolipin has only been shown to be required for pore
#   formation, prior steps in Bax activation have not been examined previously
#   (19, 20).

"""
So, is the model for Green's mode 2 inhibition at the membrane that BclXL
binds Bax with the 6A7 epitope exposed, but not inserted? Thus dissociation
from BclXL could allow Bax to return to its peripheral state, and rapidly
translocate back to the cytosol

- Every test has the observable: condition, and result

Conditions: 20nM tBid, 100nM Bax, variable liposomes:

=== Kim: (Bcl-2 Changes Topology) ====

=== Annis (Bax forms multispanning monomers) ===

Intro:

* Modulation of membrane properties is crucial to the regulation of apoptosis
  (Wolter et al, 1997).

* "For Bax to selectively target to the outer mitochondrial membrane, the
  tail-anchor sequence is required (Wolter et al, 1997; Nechushtan et al,
  1999).
* In addition, the central a5 and a6 helices (Nouraini et al, 2000) and regions
  of the amino-terminus (Goping et al, 1998) contribute to the regulation of Bax
  membrane binding.

* Typically, translocation of roughly 20% of the cellular Bax from the cytoplasm
  to the mitochondrial outer membrane is sufficient to induce apoptosis (Annis et
  al, 2001). 

* By itself, Bax can permeabilize a variety of other membranes including
liposomes (Kuwana et al, 2002) in which large oligomers of Bax have been
visualized (Epand et al, 2002). Epand RF, Martinou JC, Montessuit S, Epand RM,
Yip CM (2002) Direct evidence for membrane pore formation by the apoptotic
protein Bax. Biochem Biophys Res Commun 298: 744-749

* In contrast, based on experiments using liposomes and purified proteins, the
for- mation of lipidic pores has been proposed as yet another mechanism for
large pore formation in this case by Bax together with a large molar excess of
tBid (Kuwana et al, 2002; Terrones et al, 2004). Neither of these models is
directly pertinent to membrane permeabilization by Bax in cells in which
release of intermembrane space proteins by Bax or Bak can be induced by
100-fold substoichiometric quantities of tBid (Ruffolo and Shore, 2003).

Results

* Mitochondria from Myc-/- cells are resistant to permabilization by tBid
and Bax, even when lysate from Myc+ cells is added (Fig 1A)
* Mitochondria from Myc+ cells are sensitive to permeabilization by
tBid and Bax, even when lysate from Myc-/- cells is added (Fig 1B)

*  Indeed at 100 nM Bax, 0.4 nM tBid is sufficient to induce almost complete
release of cytochrome c from Myc+ cell mitochondria (Supplementary Figure
1). 

Read!!!
* Annis supplement
* Ruffolo and Shore
* Kuwana 2002
* Terrones 2004
* Soucie 2001
* Cartron 2003; Garcia-Saez 2004

 Good set of references for basics of pore formation

=== Lovell ===
 

* Stable, reversible tBid-Bax interaction stabilizing within ~30 minutes
* Rank order of reaction half-times
* tBid translocates to liposomes very rapidly (~10^12 molec / sec)
* Supplement: translocation of Bax to liposomes reaches ~80% within ~30 minutes


Satsoura:
* Bax not stably bound to liposomes in the absence of tBid (given Ka value)
* Saturation of Bax localization to liposomes in the presence of tBid
* Minor differences in Bax localization to liposomes of different sizes

Kale data
* Reaction half-times of different conformational changes
* Each one well-fit by a double exponential (except BH3)

Other:
* 90-95% of tBid translocates to liposomes
* "Bax jumping" whereby tBid equilibrated liposomes facilitate the permeabilization
  of no-tBid liposomes


"""
