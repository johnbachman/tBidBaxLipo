Summary/Conclusions
===================

**Multi-phasic NBD-Bax curves.**
    Justin's :doc:`data for 15 NBD-Bax mutants </datasets/nbd_bax/plate_1>`
    measured in the plate reader suggest that there are 3-4 intermediate
    states, with different residues fluorescing differently in these states.
    However, this remains to be analyzed further using data from the PTI.

**Effects of Bax-liposome titration on Bax insertion kinetics.**
    The results from Justin's :doc:`c126 titration experiment
    </datasets/nbd_bax/c126_titration>` show clearly that increasing liposome
    concentration makes Bax insertion faster, saturating at high
    concentrations.  This can be explained even with a simple solution/liposome
    :ref:`partitioning model <bax-insertion-partitioning-model>` due to the
    larger fraction of Bax peripherally bound to membranes and hence the more
    rapid production of inserted Bax. However, the simple partitioning model
    :ref:`fails to explain <bax-insertion-partitioning-model>` why increasing
    Bax concentration should :ref:`make the insertion kinetics slower
    <bax-liposome-surface>`. However, this can be addressed by :ref:`accounting
    for the saturation of the activator tBid <bax-insertion-tbid-activation>`.
    It is possible that a liposome binding site model in which Bax only
    occupies sites while it is peripherally bound and not after it is inserted,
    though this seems that it would be equivalent to the tBid case (where the
    amount of tBid is a stand-in for the "binding sites" of the liposome).

**mDivi-1 and analogs inhibit dye release in liposomes.**
    Kushnareva et al.  showed that permeabilization of outer membrane vesicles
    (OMVs) derived from mitochondria was inhibited by the addition of analogs
    of the drug mDivi-1, whose canonical mechanism is the inhibition of the
    GTPase Drp1/Dnm1, which is responsible for mitochondrial fission
    [Kushnareva2012]_. Despite the lack of detectable Drp1 in their OMVs, and
    the lack of any GTP in the system, the authors nevertheless concluded that
    the mDivi-1 analogs targeted a hypothetical "catalyst" protein responsible
    for permeabilization in OMVs. However, :ref:`experiments with synthetic
    liposomes <130830_mDivi1_analogs>` suggest that these drugs inhibit
    Bax-mediated permeabilization even in the absence of any putative catalyst
    protein. Moreover, the rank-ordering of the effect of these drugs in the
    synthetic liposome system is the same as in their OMV system, with analog H
    the most potent, then analog B, followed by the inactive analog 309s
    (actually, my experiments show that the canonical mDivi-1 is actually the
    most potent, but they did not show results for this drug in their paper).
    Interestingly, the drugs also show a base level of liposome permeabilizing
    activity, and the rank ordering of this activity is the same as their
    activity in suppressing Bax-mediated permeabilization, with analog H
    causing the most permeabilization, followed by analog B and 309s
    (interestingly, mDivi-1 itself causes the lowest basal permeabilization
    despite have the greatest effect on Bax). It is important to note, however,
    that the magnitude of the permeabilization inhibition is nearly 100% in the
    OMV system of Kushnareva et al, while only approximately 50% for mDivi-1 in
    the synthetic liposome system.  One possible interpretation of these
    results is that, first, the effect of mDivi-1 on permeabilization in both
    OMVs and synthetic liposomes is based on an interaction with either lipids,
    Bax, or tBid, and not with an unidentified catalyst protein; second, that
    the mDivi-1 analogs interact with lipid membranes to destabilize them and
    cause slow release of small-molecule marker dyes. Taken together, this
    suggests that the effect of the mDivi-1 analogs is to alter lipid
    organization or packing in such a way as to restrict or retard the
    formation of Bax pores.

**Two-phase dependence of kinetic rate on Bax concentration for heated Bax.**
    Results from :ref:`permeabilization experiments with Bax heated to 43C
    <130614_Bax_43C_titration>` show that the kinetics at each concentration
    are well-fit by a three-parameter equation, :math:`F_{max}\left(1 - e^{k_1
    (1 - e^{-k_2 t}) t} \right)`.  In this equation, the parameter :math:`k_1`
    can be interpreted as the apparent rate of permeabilization after an
    initial activation period governed by :math:`k_2`. Interestingly, plots of
    the concentration-dependence of this parameter show that it has a two-phase
    scaling with the amount of Bax: at high Bax, the rate of permeabilization
    goes up linearly with the amount of Bax, while at lower Bax concentrations
    the permeabilization rate drops off, deviating from linearity. Kushnareva
    et al. observed similar scaling in experiments with OMVs, but concluded
    that "some small deviations from linearity at low Bax concentrations
    probably reflect systematic error due to the need for long incubations
    under those conditions." [Kushnareva2012]_ In my experiments this deviation
    is highly reproducible and is robust to the choice of model used. To
    explain the two-phase nature of this kinetic scaling, I considered
    oligomerization, transient binding to lipid membranes, auto-activation of
    Bax, and other mechanisms, but found in systematic modeling studies that
    none of these were able to reproduce this feature of the data. However,
    when I allowed for two lipid binding sites, one saturable and one
    non-saturable, the fit to the data greatly improved. In this hypothesis the
    measured rate :math:`k_1` is proportional to the amount of transiently
    membrane-bound Bax, :math:`Bax_m`. Non-saturable binding predicts that
    :math:`Bax_m` will increase linearly with total Bax; saturable binding
    predicts that at high Bax, :math:`Bax_m` will plateau. A two-site model
    predicts that there is a high-affinity, saturable site and a low-affinity,
    non-saturable site. This type of model yields steady-state concentrations
    of :math:`Bax_m`, and dye release rates, that have the observed two-phase
    characteristics.

**Deterministic approximations of dye release break down for some mechanisms.**
    See the :ref:`introduction to the section on stochastic permeabilization
    models <stochastic_models_intro>`.

**tBid titration shows no inhibition at high tBid.**
    Need to add figures and data here.

**Biphasic character of dye release/pore curves.**
    This was most visible in the :doc:`first round of dye release curves
    </datasets/ants/tbid_bax_lipo_1>`. However this effect, though still
    visible, was much reduced in the :doc:`second round of curves
    </datasets/ants/tbid_bax_lipo_2>`. Hypotheses considered included tBid-Bax
    inhibition or Bax depletion and subsequent recycling through dissociation.
    Another explanations to consider could be tBid dissociation (jumping).

**Biphasic character of NBD-Bax insertion curves (by PTI).**
    This was most visible in the curves that Justin measured in the PTI,
    especially the c62 mutant. A possible explanation for the c62 signal
    was tBid-Bax inhibition coupled with a fluorescence change that depends
    on a protein-protein interface between tBid and Bax (and/or between
    Bax and Bax).

