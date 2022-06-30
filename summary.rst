####################
Overview and summary
####################

This part of the best practice guide provides an overview and summary
of best practices regarding QM/MM simulation of biomolecular systems
considered to be important by the organisers and speakers of the
`BioExcel Virtual Workshop on Best Practices in QM/MM Simulation of
Biomolecular Systems
<https://bioexcel.eu/events/virtual-workshop-best-practices-in-qm-mm-simulation-of-biomolecular-systems/>`_.

===============
1. Introduction
===============

QM/MM simulation is among the most versatile tools to model systems,
because it relies on an *ab initio* description of the most important
part of the system. By performing such simulations one has access to
properties that sometimes escape experimental probes, like transition
states, short-lived intermediates, conformational itineraries,
electronic states and so on. QM/MM has been applied to investigate:

* Reaction mechanisms, including rate selectivity and the effect of
  point mutations
* Ligand protein interactions and binding affinity estimations
* Crystal structure refinement, in particular for systems with exotic
  co-factors

However, performing a QM/MM study is not as straightforward as for
other computational approaches employed to investigate biological
systems and care has to be taken when setting up a QM/MM simulation.


================
2. Prerequisites
================

QM/MM methods are computationally very demanding compared to other
approaches and for this reason they have traditionally rarely been
considered when describing biological processes. However, there are
cases where their use can be very valuable. For example, by employing
QM/MM methods at a suitable QM level of theory, one can nowadays make
reliable predictions about activation barriers and reaction mechanisms
in enzymes. There are a lot of mechanisms in biochemistry textbooks
that are probably wrong because those mechanisms were written down at
a time when all one could do was to speculate about them on the basis
of kinetics and perhaps mutagenesis experiments. Since enzymes often
use some quite unusual chemistry, intuition and well-informed guesses
alone can unfortunately lead to inaccurate conclusions. A suitable,
well-executed QM/MM calculation is therefore advantageous as it can
provide a sensible mechanistic prediction of an enzymatic mechanism.

As probably in any computational study of biological systems, the
first preliminary step before setting up any QM/MM simulation is to
study the literature about the system under investigation and the
properties you are going to measure. For example, the apparent
activation energies measured so far over all known enzymatic classes
range between 5 to 25 kcal/mol, and more than 80% of energies are
between 14 and 20 kcal/mol. Therefore, if you obtain a value outside
those intervals, you should suspect that something was wrong with your
calculations.


======================================
3. Choosing method and level of theory
======================================

Before setting up a QM/MM simulation, it is important to remember that
QM/MM is not a “black box”. Care must be taken in choosing the
suitable QM/MM scheme and the suitable level of theory, in particular
for the QM part. Nowadays, especially when large QM regions (order of
hundreds atoms) are involved, the density functional theory (DFT) is
by far the most common level of theory employed to treat the quantum
part, because it usually provides a good tradeoff between accuracy and
computational cost.

A downside of DFT is that it is not systematically improvable, i.e. we
do not know *a priori* if the result with one functional should be
better or worse than the results with a different functional. When
this is a requirement to obtain a mechanistically sound prediction, then
one must go beyond DFT and employ some post-Hartree Fock
approach. Among the post-Hartree Fock approaches applied to large
systems, the family of coupled cluster methods (e.g. CCSD(T)) are
still too computationally expensive for typical biological
applications and they have been rarely employed so far, in spite of
their state-of-the-art accuracy. A better balance between high
accuracy and computational efficiency is provided by approaches based
on the second order Møller–Plesset perturbation theory (MP2), such as
the spin-component scaled MP2 (SCS-MP2) as proposed by Grimme
(`J. Chem. Phys. 118, 9095 (2003)
<https://doi.org/10.1063/1.1569242>`_), in which the accuracy of MP2
calculations is significantly improved by semi-empirically scaling
the opposite-spin and same-spin correlation components with separate
scaling factors.

In general, the selection of the method to employ, i.e. the
Hamiltonian which is going to describe your system, will require you
to carefully consider:

* if the Hamiltonian correctly describes all the energy contributions
  involved in the phenomenon you want to investigate;
* the availability of software that implements it;
* the computational resources at your disposal.

Within DFT, the choices of functional and basis set are crucial to
obtain a correct estimation of the physical quantity to be predicted.
The basis set should contain a sufficient number of elements to have a
correct description of the wave function, but at the same time it
should be as small as possible in order not to make the calculation
unfeasible. This can (and should!) be validated via convergence tests
on some relevant physical quantity, such as the total energy. In
general, structural properties require smaller basis sets than
properties based on the estimation of energies.  In cases where one is
performing a calculation using a localised basis set, a common
recommendation is to always employ polarisation functions (the minimum
requirement is a valence double-zeta basis set). Diffuse functions are
employed less frequently in QM/MM schemes due to their large extension
and because of possible artifacts originating from interaction of the
electronic density with partial charges in the MM region. If diffuse
functions are needed, one can employ only diffuse functions associated
with the atoms in the central area of the QM region, far from the
QM/MM boundary.

The selection of a suitable functional on the other hand is not so
straightforward. A thorough review of the literature about the system
under investigation could help make an educated guess regarding the
best functional to employ. Even if your particular system has not been
investigated yet, it is highly probable that the chemistry involved in
it has been widely studied with DFT, and that several benchmarks exist
in the literature that could help make an appropriate guess. However,
in general, especially when the information from the literature is
scarce, an *a posteriori* validation of the functional on similar
properties and/or similar (smaller) systems is recommended (see
section 5) before embarking on expensive QM/MM calculations. The
validation can be done with respect to experimental data, if
available, or by modelling a minimal system that contains at least the
relevant part (e.g. the active site) and performing on it a high-level
calculation (e.g. CCSD(T)/CBS) to be used as a reference for the
energies. Moreover, in the case of enzymatic reactions it is in
general a very good test to compare reaction energetics between gas
phase, water solvation only, and the complete enzyme environment: if
the level of theory is adequate, one should see a clear catalytic
effect in the protein environment.

Having mentioned solvation, as a general remark explicit solvation
usually provides better agreement with experimental results compared
to any implicit or continuum solvation models. However it also brings
higher computational costs.

A popular and computationally inexpensive strategy that often allows
improving both structure and energy predictions from “pure” DFT
approaches is to add the so-called dispersion corrections. These are
suitable additional terms added to the functionals that can yield
results in much better agreement with experimental data. In fact, in
biomolecules and similar systems dispersion competes significantly
with other effects, while DFT is known to provide an incomplete
treatment of dispersion, which can adversely affect its accuracy.

Once the QM level of theory and the MM description have been decided,
one should focus on the selection of the most suitable QM/MM scheme to
employ. A difficulty with QM/MM approaches is that there exist so many
variants and that the details of these are seldom discussed. For
example, QM/MM approaches can use either subtractive or additive
schemes (`Senn and Thiel, 2009
<https://doi.org/10.1002/anie.200802019>`_).

In a subtractive scheme, three separate calculations are performed:

#. QM calculation with the QM region
#. MM calculation with the entire system 
#. MM calculation with the QM region

Then the total energy of the QM/MM system is estimated as the sum of
the first two terms subtracted from the third term in order to avoid
counting twice the interactions within the QM subsystem. The advantage
with this approach is its simplicity: it automatically ensures that no
interactions are double counted and it can be set up for any QM and MM
software (provided they can write out energies and forces), without
the need of any modification of the code. The typical example of a
subtractive scheme is ONIOM (`Svensson et al., 1996
<https://doi.org/10.1021/jp962071j>`_) or ComQum (`Ryde, 1996
<https://doi.org/10.1007/BF00402823>`_).

In an additive scheme, there are no calculations at different
resolutions taking place separately. Instead, the potential energy for
the whole QM/MM system is a sum of three contributions:

* QM energy terms;
* MM energy terms;
* QM/MM coupling terms.

In contrast to the subtractive schemes, the interactions between the
particles in the QM region and the classical atoms in the MM region
are treated explicitly in through QM/MM coupling terms. In this case,
it is up to the developer to ensure that no interactions are omitted
or double counted. Therefore, an additive scheme requires special MM
software, in which the user or developer can select which MM terms to
include. The advantage of the additive QM/MM scheme is that no MM
parameters for the QM atoms are needed, because those energy terms are
calculated by QM. For this reason, additive schemes are gaining more
and more popularity among the QM/MM user community and nowadays are
the preferred QM/MM approach in biological applications.


The interaction between the QM and MM regions is typically dominated
by electrostatics. This interaction can also be considered at
different levels of approximation. These approaches can be classified
as either mechanical embedding, electrostatic embedding or polarized
embedding. In mechanical embedding - the simplest electrostatic
coupling scheme - the QM-MM interaction is calculated at the MM level
and the electronic wave function is computed for an isolated QM
subsystem. In other words, the MM environment cannot induce
polarization on the electron density in the quantum region. Moreover,
since the QM region is usually the site of the reaction it is likely
that during the course of the reaction the charge distribution will
change, resulting in a high level of error if a single set of MM
electrostatic parameters is used to describe it. For those reasons,
mechanical embedding is nowadays no longer recommended for modeling
reactions in biochemical macromolecules. Instead, the large majority
of the state-of-the-art QM/MM codes employ electrostatic embedding, in
which the electrostatic interactions between the QM and the MM
subsystem are treated at QM level and handled during the computation
of the electronic wave function. This is done by including the atomic
partial MM charges in the Hamiltonian used to solve the quantum
problem, i.e. the Hamiltonian depends on both the classical partial
charges and the quantum charge density.

Increasing further the level of sophistication means including in the
model also the polarizability of the MM atoms. This is done through a
polarization embedding scheme, in which both regions - QM and MM - can
mutually polarize each other. Although this last embedding offers the
most realistic electrostatic coupling between the quantum and the
classical regions, polarizable force fields for biomolecular
simulations are not so effective yet. Therefore, despite progress in
the development of such force fields, QM/MM studies with polarizable
MM regions are so far not so popular.

Moreover, in several systems it has been shown that QM/MM approaches
employing electrostatic embedding allow reaching almost the same
accuracy as obtained by full QM treatments of the entire system
(computationally much more expensive). However, to reach such an
accuracy, a user usually needs a rather expert knowledge of the
methods and their approximations, and significant experience in
defining suitable domains/regimes within which to apply them. It is
worth mentioning here that recently schemes have been developed to
make it easier and more straightforward to access high-level,
potentially chemically accurate QM/MM treatments. One example is the
projector-based embedding scheme developed by Manby, Miller and
co-workers (`JCTC 8(8), 2564-2568 (2012)
<https://doi.org/10.1021/ct300544e>`_). Projector-based embedding
allows for rigorous high-level treatment of a small region (e.g. with
coupled cluster or MP2 level of description) within a larger DFT
region, and can be incorporated within a QM/MM framework. As mentioned
before, the selection of a suitable functional within a DFT approach
is crucial: reaction energies and barriers can significantly depend on
the employed DFT functional. However, such a dependence disappears
within approaches like the projector-based embedding, and energies and
the barriers can be then predicted with a high degree of confidence.

Most of the electrostatic embedding QM/MM approaches introduce long
range electrostatic cutoffs of some kind in their schemes in order to
reduce the computational load. Therefore, one should have an idea of
the extent to which long-range interactions are important to correctly
describe the phenomenon under investigation in one's system of
interest. For example, in enzymes electrostatic long-range
interactions typically cover a range from 7 to 20 Å and therefore you
have to be careful to include them in the QM/MM modelling when
planning simulations.



==================
4. QM/MM modelling
==================

The usual starting point in preparing a QM/MM simulation is a PDB
structure coming either from an X-ray, NMR or cryo-em experiment, or
as the output of a homology modelling or other theoretical protein
structure prediction approach. It is not granted that such an initial
structure is immediately usable. The initial PDB could have missing or
incorrectly identified atoms/residues, and the protonation state of
histidines and/or ligands could be incompatible with their chemical
environment, to name just a few potential issues. The initial
structure must be critically evaluated, and all such issues need to be
addressed before proceeding. For example, to assess the protonation
state of ionizable groups in protein and protein ligand complexes we
need in principle to find their pKa’s. If already known, pKa values
can be found in the PubChem database
(`https://pubchem.ncbi.nlm.nih.gov/
<https://pubchem.ncbi.nlm.nih.gov>`_).  Alternatively, theoretical
methods have been developed to make reliable estimations, and some of
them have been implemented in servers and tools, such as MarvinSketch
(`https://www.chem.uwec.edu/marvin/
<https://www.chem.uwec.edu/marvin>`_) and PROPKA
(`https://github.com/jensengroup/propka
<https://github.com/jensengroup/propka>`_), the latter being
particularly popular among computational enzymologists.

After obtaining the correct initial structure, one needs to think
about how to partition the system, i.e. what the QM region should
contain. Of course, the QM part needs to include all the atoms
involved directly or indirectly in the chemical process under
investigation. Monitoring the possible interactions between individual
residues occurring during a preliminary force field-based simulation
using classical Molecular Dynamics may help with this step. If time
and computational resources allow it, one should test QM regions of
increasing size. It is a well known problem that the energies obtained
from both QM cluster (i.e. full QM description without any
partitioning) and QM/MM calculations depend on the size of the QM
part, and that these energies converge slowly with respect to the size
of the QM region. Even worse, they do not necessarily converge in all
cases: the impact on accuracy of some of the underlying limitations of
the QM/MM approach will in fact increase in severity with increasing
system size. For example, as we scale up the size of the QM region
within a QM/MM calculation, the QM/MM boundary becomes larger too and
if there is excessive polarization of that boundary, that error will
increase as well. Therefore, simply increasing the size of the QM
region is not necessarily a good test of a QM/MM calculation and that
is especially the case if you are dealing with a heterogeneous
environment like an enzyme, where the charge of the system may differ
as you increase the size.  Moreover, quite often the variation of the
energy values is less dependent on the size of the QM region than it
is on the choice of the DFT functional, as has been shown for example
in several reaction energy and barrier calculations. Therefore, one
should probably pay more attention to the suitable QM treatment in
describing the system than the size of the QM region.

In cases where the convergence with QM region size has been studied
and successfully obtained, typically something in between 500 and 1000
atoms were needed to obtain convergence to reasonably accurate
results! Dealing with such large QM regions is a challenging task for
most quantum codes and requires significant computational
resources. Luckily, QM/MM energies (and also structures) normally
converge faster than the QM cluster ones. General recommendations to
get stable energies (and structures) in a protein system is to include
in the QM region:
- Neutral groups/residues up to 4-5 Å from the active site
- The largest number (ideally all) of the charged groups that are not
on the surface of the protein, i.e. the buried ones.

In addition, the “junction” atoms, i.e. the atoms at the border
between the QM and the MM boundary should be kept as far as possible
from the active site.

In fact, when the boundary between the QM and MM regions cuts a
covalent bond connecting a quantum atom to a classical atom, care has
to be taken in solving the quantum problem, i.e. in calculating the
wave function associated with the QM region. Biological systems have a
large content of sp3-hybridized C-C bonds and due to their symmetry
and lack of polarization, those bonds are the best choice of locations
where to place the QM/MM boundary. However, a straightforward cut
through the QM/MM boundary would create one or more unpaired electrons
in the quantum subsystem. In reality, these electrons are paired in
bonding orbitals, with electrons belonging to the atom on the MM
side. However, now those electrons do not exist in the MM region, due
to the artificial partitioning and the lower level of resolution we
decided to use in this region. In the literature a number of
approaches have been proposed to remedy the artefact that originates
from such open valencies, which can generally be classified into three
categories:

#. The link-atom approach, in which additional link-atoms, generally
   hydrogen atoms, are introduced at an appropriate position along the
   bond vector to saturate the dangling valences of the quantum
   region;
#. Alternatively, it is possible to use the link atom pseudopotential
   approach, which consists in introducing in the QM system a
   description of the classical atom at the border bonded to the
   quantum atom, through a special pseudopotential with the required
   valence charge.
#. The last category involves the use of bonding hybrid orbitals,
   including the hybrid orbital method, the local self-consistent
   field method, the generalized hybrid orbital method, and the frozen
   orbital method.

All the above-mentioned methods are more or less equivalent if
carefully applied according to the specifications and recommendations
of the developers. For example, the link atom pseudopotentials require
constraining the bond distance appropriately. Therefore, the choice of
which category and also which specific approach to use often depends
on the particular implementation available in the software we are
planning to use for QM/MM calculations.

=========================================
5. Simulation protocol and QM/MM workflow
=========================================

QM/MM calculations that involve only single conformations are
computationally less demanding and can be particularly suitable to
identify *structural* features of reactive conformations. For example,
single conformation QM/MM calculations may can efficiently enable
better understanding of structure-activity relationships if
computational resources are limited. As a general remark, QM/MM single
conformation approaches are a good way to calculate structures in a
protein system, especially when combined with experimental data, while
properties based on energies are much harder to obtain since proper
sampling and solvation is required. In fact, the subset of
catalytically competent conformations can be significantly small in
comparison with the full conformational landscape, and X-ray
structures are usually a good starting point for reactivity
study. However, a re-optimization of the initial structure(s) is
usually recommended.

Since geometry optimization, which allows the extraction of structural
information, is typically computationally more demanding than wave
function minimization, from which one can obtain energetic
information, it is common to use a lower level of theory for the
former and a higher level for the latter. Although this strategy can
usually be employed with confidence, we should always bear in mind not
to mix very different Hamiltonians in performing these two sequential
steps, for example combining very different DFT functionals, e.g. a
GGA family functional like PBE and a hybrid functional like B3LYP. At
the very least one should double check any conclusions arising from
such “dangerous” mixing by testing different combinations.  In
general, techniques that allow more efficient exploration of the
conformational space of an enzyme during catalysis can provide a more
dynamic picture of the PESs associated with the reaction. For this
reason, it can sometimes be beneficial to perform preliminary
classical (i.e. force field-based) molecular dynamics simulations of
the system, select representative snapshots according to predetermined
adequacy criteria, and to then apply the QM/MM single-point protocol
to these snapshots.

If in contrast temperature and entropic effects are important in the
process under investigation and cannot be neglected, then sampling
through a QM/MM molecular dynamics approach is required. Even then
however, the system is usually still thermally equilibrated at
classical level before performing a QM/MM molecular dynamics
simulation. Classical molecular dynamics allows sampling of longer
time scales not accessible within a QM/MM molecular dynamics
simulation. Therefore, this equilibration step is particularly
relevant for those systems which are very flexible, such as complexes
formed by carbohydrate-active enzymes, and for those systems whose
initial structure can not represent a configuration with sufficiently
high probability to be found in the corresponding ensemble.

On moving from the classical to the QM/MM description, another thermal
equilibration is required before starting the QM/MM production
run. This is because when the level of theory is changed, the forces
felt by the atoms change as well. This second equilibration step
allows the system to relax and adapt to the new QM/MM Hamiltonian.

Enhanced sampling methods can be used when standard classical or QM/MM
molecular dynamics approaches do not manage to sample enough
configurations or parts of phase space in a reasonable time due to the
presence of large free energy barriers or if doing so would be much
beyond the computational resources at one's disposal, which may the
case especially given how large systems biological systems often
are. Most enhanced sampling methods fall into two categories:
collective variable-based and collective variable-free methods.

A direct and effective idea to accelerate the thermodynamics
calculation is to modify the potential energy surface by adding a bias
potential to the Hamiltonian of the system, thereby decreasing the
energy barrier hence increasing the sampling of transition
regions. This is the strategy followed by methods such as the original
and widely used umbrella sampling (`Torrie and Valleau, 1977
<https://doi.org/10.1016/0021-9991(77)90121-8>`_) and the popular
metadynamics (`Laio and Parrinello, 2002
<https://doi.org/10.1073/pnas.202427399>`_).

The second category include methods such as parallel tempering
(`Swendsen and Wang, 1986
<https://doi.org/10.1103/PhysRevLett.57.2607>`_) and replica exchange
molecular dynamics (REMD, `Sugita and Y. Okamoto, 1999
<https://doi.org/10.1016/S0009-2614(99)01123-9>`_), which alter the
canonical probability distribution to a distribution that induces a
broader sampling of the potential energy.

The collective variable-based methods need to use predefined reaction
coordinates or collective variables (CVs), i.e. low-dimensional
functions of electronic and atomic degrees of freedom of the system
being simulated, which should refer to the variable describing the
slow motion in the process of interest. Appropriately selecting the
CVs is crucial for such methods to correctly work. For example, since
the transition states are not identified correctly by the CVs, forward
and backward transitions might follow different paths, leading to
hysteresis in the estimated free energy, or the bias potential could
show large fluctuations, and the free energy estimate does not
converge. However, it is well known that the proper reaction
coordinates are not easily identified for many systems. In describing
a reaction pathway, a CV should include all the bonds which form or
break during the process. Low-energy vibrations and variables in most
cases do not significantly affect free energy landscapes, therefore
they should not be involved in the CVs. Literature research can be
very beneficial for this task because a lot of different CVs have
already been tested and used to describe a plethora of different
phenomena. When possible, a trial and error approach using simpler
model systems or suitable enhanced sampling parameters to quickly get
a rough, low-resolution estimate of the underlying free energy
landscape (e.g. choosing hills with large heights in metadynamics) can
be useful in identifying the proper reaction coordinate(s).

================================================================
6. Analysis, interpretation and validation of simulation results
================================================================

After performing QM/MM modelling one should check the validity of
results obtained before publishing. There are three major things that
need to be established and validated:

#. QM theory level is sufficient and predictive for your system
#. Structures (models) and energies of reactants, products and transition states
#. Reaction pathways should represent mechanism and be predictive

The level of QM treatment of the system should be sufficient for
predictivity but typically it can not be too high because of limited
time and computational resources. At very high levels it is usually
possible to calculate only energies (enthalpies) for a handful of
structures. A systematic way to check validity in that case is to
obtain energies at a lower level of theory and characterise them with
higher-level methods. It is also good to consider structures obtained
with free energy methods and characterize them with static methods
(i.e. geometry optimizations) at a higher level of theory.

As mentioned in previous parts, models can be divided into static
(obtained with geometry minimizations) and dynamic (obtained with
MD). Static methods can, in the case of enzymes, give only good
estimation of enthalpies, so one should validate whether entropy has
only a minor effect on reactivity. That could be estimated, for
instance, with QM/MM dynamics using a lower level of theory. Note as
well that the effect of temperature on the reaction can be very
effectively estimated with QM/MD methods as well, e.g. using
metadynamics or umbrella sampling.

Substrate conformations can be wrong due to forcefield inaccuracies if
you are using starting structures extracted from classical MD
simulation. In such cases it is valuable to use QM/MD to relax the
starting geometry. Another thing to validate in the model is the
reaction pathway and variables used to represent it in the
simulation. It is always good to start from a simple model of your
system in solution or cluster (with polarizable continuum) and test
various reaction pathways beforehand. In addition, good collective
variables that represent the reaction typically include all breaking
and forming bonds. If you are using either position restraints or
constraints in your model, remember to check their effect on the
reaction energetics as they may affect final results heavily and
should not be used for production of the final reaction profile.

Reaction pathway validation can be a major obstacle. First of all,
reaction space in enzymatic reactions can be very complex, with many
individual steps. Second, sampling time and choice of initial
structures are crucial as in many cases the reactive conformation of
the enzyme-substrate complex is not the most populated during
sampling. Check as many initial structures as possible. Sometimes the
Michaelis complex in the enzymes could have a distorted geometry that
is not energetically favorable, but resembles a transition state. In
dynamical methods sampling time is very important, check that your
profile is converged by increasing sampling time and confirming that
further simulation does not change your free energy profile. In case
of very fast enzymes, the Michaelis complex formation constant plays a
key role. For instance, if one is trying to model the effect of a 
mutation on enzyme activity then it should be checked how the mutation
affects not only reaction barriers but also binding free energy of the
substrate. One way to do that is to compare reaction energetics
between solution (cluster model) and enzyme (QM/MM model).

Final results require validation as well. This is typically done using
existing experimental data, if available, and/or by making predictions
that could be verified against experiment. In the case of well-known
systems, validation against experimental data should be done on the
whole array of available results. In other words, results obtained
through simulation should be in line with all available experimental
data. Do not cherry pick only those experiments which correlate well
with your results, omitting others that do not align with your
hypothesis from simulation. Do not forget that scientific results
should be falsifiable, meaning that there should always be some
prediction as a result of your work, which could be proven (or
rejected) with an experiment. Some effects do not contribute directly
to enzyme catalysis. For example, kinetic isotope effect (KIE) may not
be included directly in the energy of the transition state, however
sometimes its effect (i.e. in case of proton tunnelling) should be
properly included into calculation of the catalytic constant in order
to obtain a predictive model. If mutation experiments are used for
validation one should be aware that they could affect structure and
energy of the enzyme-substrate complex, thus changing not just the
barrier but also the Michaelis constant (for the binding pose and/or
equilibrium). The same validation should be done for inhibitors and
non-standard substrates studied. In the case of a static model, to
check that it did not fall into the local minima, try to build the
profile both from reagents to products and back from products to
reagents, if your model is correct the profiles should be identical.


  
