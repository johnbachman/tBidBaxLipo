from pysb import *
from pysb.macros import reversibly_translocate

# Common to all models
Model()

Monomer('Bax', ['c3', 'c62', 'c120', 'c122', 'c126', 'c184'],
    {'c3': ['s', 'm'],
     'c62': ['s', 'm'],
     'c120': ['s', 'm'],
     'c122': ['s', 'm'],
     'c126': ['s', 'm'],
     'c184': ['s', 'm']})

Observe('Baxc3', Bax(c3='m'))
Observe('Baxc62', Bax(c62='m'))
Observe('Baxc120', Bax(c120='m'))
Observe('Baxc122', Bax(c122='m'))
Observe('Baxc126', Bax(c126='m'))
Observe('Baxc184', Bax(c184='m'))

# Takes a monomer pattern as an argument
def all_parallel():
    for site in Bax.sites:
        reversibly_translocate(Bax(), 's', 'm', locname=site)

def linear_pathway(num_phases):
    """ Partition the sites into n groups.
    Rule 1: from all sites in solvent to insertion of sites in first group.
    Rule 2: from all sites in first group inserted to this plus all sites in second group.
    Rule 3: from all sites in first and second groups inserted to this plus all sites in third group
    Rule 4: ...etc."""

    # Should return list of (n) lists of sites (strings) 
    group_list = partition(Bax.sites, n)
    
    # Iterate over site groups
    all_prev = {}
    for i, group in enumerate(group_list):
        # Make the rule
        # Make the left hand side dict {'c3':'s', 'c62':'m'...}
        # Iterate over all previous groups
        for prev_group in enumerate(group_list[0:i]):
            # For each group, iterate over sites
            lhs_sites = [site for site in prev_group]
              
                
        Rule('r' + str(i),
             Bax

    """
    For each phase, choo
    reversibly_translocate(Bax(), 's', 'm', locname='c62')

    Rule('b1', Bax(c3='s', c62='s', c120='s', c122='s', c126='s', c184='s') <>
               Bax(c3='s', c62='s', c120='s', c122='s', c126='s', c184='s'), kf, kr)

 ) <> Bax(c62='m', c3='m'), kf, kr)
    Rule('b2', Bax(c3='m', c120='s') <> Bax(c3='m', c120='m'), kf, kr)
    Rule('b3', Bax(c62='u', c3='s') <> Bax(c62='m', c3='m'), kf, kr)
    Rule('b4', Bax(c62='m', c3='s') <> Bax(c62='m', c3='m'), kf, kr)
    Rule('b5', Bax(c62='m', c3='s') <> Bax(c62='m', c3='m'), kf, kr)
    Rule('b6', Bax(c62='m', c3='s') <> Bax(c62='m', c3='m'), kf, kr)
    """



"""
Way to approach this:

- A protein has a discrete set of conformational/insertional states.

- Each state has associated with it a set of fluorescence magnitudes for each of the
  NBD positions.

- The signal coming from any given NBD mutant experiment will therefore come from
  a combination of both the amount (concentration) of protein a state and the fluorescence
  magnitude associated with that state. Specifically, the total NBD fluorescence
  for mutant i will be the dot product of the concentrations of each state
  with the NBD intensities associated with each state.

- A "pathway" will therefore consist of a set of transitions between these states; the
  transitions are defined by rate parameters.

For example, a two-state model:

S1 (intensity vec) <-- kf/kr --> S2 (intensity vec)

For this simple example one could imagine that the intensities of all residues is 0 in S1,
and some value in S2. As the system reached equilibrium between S1 and S2, the total
fluorescence would be determined by the ratio of the concentrations of S1 and S2.

In fact, in this specific case the dynamics would be of the form
(1 - e^-(kf+kr)t)* (kf/kf+kr) S1_0

For a three-state model:

S1(0 vec) <--> S2(vec1) <--> S3(vec2)

- However, this leads to a large number of free parameters. 5 per state, with two rate
parameters for each transition. Then again, with a double exponential fit, you have two
fmax() values for each mutant, plus two k()s for each as well, yielding a greater
number of parameters!!

In general, using a state model, for n states (not counting the zero state), you would get
5*n intensity constants + 2*n = 7n parameters to be fit.

One way to constrain the number of free parameters would be to assume that each residue
is one or two-phased; it could then have intensities assigned to at most one or two states
(all subsequent states would have the same intensities as the prior state unless allowed
to have a new intensity). 

Consider first a one (non-zero) state, one-phase model. That is, this would just consist
of S1 <-> S2, and only S2 wold have non-zero intensities associatd with it. This would
lead to single-exponential kinetics of each residue.

Now consider the two state model, but keeping it one-phase. That is, each intensity can
be changed only once along the reaction pathway. So all would start out non-zero. The N
(for NBD) residues could be distributed among the k states. The problem is basically one 
of how to partition the N things into k states. Basically, I have N! ways of arranging the
parameters, and can put the k-1 partition(s) in one of N places. But once partitioned,
I am overcounting, since the order in each bin doesn't matter.

So for nothing in the first bin, all 5 in the second, I have to divide by 5!. For one
thing in the first bin, four in the second, I have to divide by 1*4!. For the next,
I divide by 2!3!. Something is wrong here.



Measurements:

3C, 62C, 120C, 122C, 126C, a9

One model:

Everything happens in parallel:
3C -> 3Ci
62C -> 62Ci
... etc.
Rule for each site flipping independently.
If this were the case, you would get a whole bunch of guys at equilibrium
with mismatched localizations though!

Everything happens simultaneously (only two-states)
S0 (all solvent)
S0 -> S1 (3C, 62C, 120C, 122C, 126C)
You would do this by grouping the sites together in a single rule.

There are multiple states:
S0 -> S1 (3C, 62C) -> S2 (120C, 122C) -> ...
In this case the "topology" is the number of states and assignment of measurements to states
(how many combinations?)
Number of steps: 2-6
For a given number of steps n
can divide 6 things among n bins, i.e. 6 choose 2

Parallel processes
     S1 (3C, 62C, etc.)
     /
S0 <
     \   
     S2 (3C, 62C, etc.)

Two branches can have different lengths

Branched
"""


def separate_species():
    Monomer('c3', ['loc'], {'loc': ['s', 'm']})
    Monomer('c62', ['loc'], {'loc': ['s', 'm']})
    Monomer('c120', ['loc'], {'loc': ['s', 'm']})
    Monomer('c122', ['loc'], {'loc': ['s', 'm']})
    Monomer('c126', ['loc'], {'loc': ['s', 'm']})
    Monomer('c184', ['loc'], {'loc': ['s', 'm']})

    Initial(c3(loc='s'), Parameter('c3_max', 5))
    Initial(c62(loc='s'), Parameter('c62_max', 5))
    Initial(c120(loc='s'), Parameter('c120_max', 5))
    Initial(c122(loc='s'), Parameter('c122_max', 5))
    Initial(c126(loc='s'), Parameter('c126_max', 5))
    Initial(c184(loc='s'), Parameter('c184_max', 5))

    Observe('c3', c3())
    Observe('c62', c62())
    Observe('c120', c120())
    Observe('c122', c122())
    Observe('c126', c126())
    Observe('c184', c184())

