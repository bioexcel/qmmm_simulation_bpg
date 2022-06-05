####################################
Detailed advice on individual topics
####################################

This part of the best practice guide lists individual topics - all
important issues regarding QM/MM simulation of biomolecular systems -
that arose during the `BioExcel Virtual Workshop on Best Practices in
QM/MM Simulation of Biomolecular Systems
<https://bioexcel.eu/events/virtual-workshop-best-practices-in-qm-mm-simulation-of-biomolecular-systems/>`_

For each topic a time-indexed link is provided into a recording on
YouTube of the workshop presentation in which the speaker in question
shares and discusses best practice regarding that topic. In addition
to the individual speaker presentations linked to below, many topics were
also covered in the workshop's `Panel Discussion
<https://www.youtube.com/watch?v=iF05I-r6YW8&list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso>`_.


===================================================
How can QM/MM best be used in biomolecular research
===================================================

----------------------------------
When to use QM/MM, and when not to
----------------------------------

* **Enzymatic reactions, transition states inside proteins** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=172>`_)


* **Enzymes may use quite unusual chemistry, well-educated guesses on the mechanisms are often wrong, which simulation can reveal** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1258>`_)


* **To predict enzyme selectivity** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2193>`_)


* **QM/MM is very good for modelling reactivity of covalent bond breaking/formation** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2873>`_)


* **Ligand-protein interactions and binding affinity calculations** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3547>`_)


* **If a reaction involves proton transfer, consider using path integration methods to account for tunnelling** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1247>`_)


* **Model effects of point mutations** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2188>`_)


* **Do QM/MM modelling when access to experimental data is difficult, e.g. at atomic resolution in transition states, short-lived intermediates, etc.** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=431>`_)


* **QM/MM simulation can be used for crystal structure refinement** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2997>`_)


---------
Key steps
---------

What are key steps when performing QM/MM simulation?

* **Typical steps in a QM/MM simulation protocol** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1841>`_) **:**
    * **Obtain PDB structure (e.g. of an enzyme complex)**
    * **System preparation**
    * **Thermal equilibration (classical MD)**
    * **Choose QM and QM-MM coupling treatments (region size, level of theory, etc.)**
    * **Thermal re-equilibration (QM/MM MD)**
    * **Choose collective variables**
    * **Reaction simulation (static or dynamic, e.g. metadynamics)**
    * **Analyze the free energy landscape**
    * **Determine reaction mechanism** 

----------------------------------------
Static modelling vs dynamical simulation
----------------------------------------

When to use static modelling (single point energies, nudged elastic band, etc.) vs dynamical simulation (umbrella sampling, metadynamics, etc.)


* **Sampling is required, however sampling time to achieve meaningful results can be greatly reduced** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3021>`_)


* **Free energy simulations are more reliable for enzymes** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2344>`_)


* **The subset of catalytically competent conformations can be significantly small in comparison with the full conformational landscape, and x-ray structure is a good starting point for reactivity study. For example, single conformation QM/MM calculations may better enable understanding structure-activity relationship if computational resources are limited** (`listen to Maria João Ramos on YouTube <https://youtu.be/XIHMcR_tR7E?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1839>`_)


* **Static techniques that allow efficient exploration of the conformational space of the enzyme during catalysis can provide insights into the dynamical picture through the PESs associated with the reaction**:
    * **Single conformation calculations** (`listen to Maria João Ramos on YouTube <https://youtu.be/XIHMcR_tR7E?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1802>`_)
    * **Averaging over an ensemble of conformations** (`listen to Maria João Ramos on YouTube <https://youtu.be/XIHMcR_tR7E?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2481>`_)

|

============================================
Structure / model preparation and validation
============================================

------------------------------------
Structure preparation and validation
------------------------------------

How to set up and validate your model structure (X-ray, NMR,
cryo-em, homology modelling, protein structure prediction), including
how to consider missing residues, atoms, protons (!), rotameric
and tautomer (again protons!) states, (missing) waters, and substrate. 


* **It is crucial to consider conformational behaviour of the substrate inside the enzyme pocket. Sample with classical or QM/MM MD.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2371>`_)


* **Evaluate critically your initial PDB structure. If needed, re-optimise the crystal structure.** (`listen to Ulf Ryde on YouTube <https://youtu.be/aQdjC-W9Wy4?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=302>`_)


* **Conformation of substrates in classical MD can be wrong due to forcefield parameterization. QM/MM can be used for sampling to clarify conformation stability.** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=795>`_)

* **It is good to check if your simulations can reproduce features from the crystal structures, like distortions of the substrate or amino acids. Simulations can be used for further crystal structure refinement** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2997>`_)


* **Suggestions to overcome the local minima problem in the case of proteins** (`listen to Ulf Ryde on YouTube <https://youtu.be/aQdjC-W9Wy4?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2009>`_) **:**
    * **Run back and forth from starting state to final state until convergence**
    * **Optimise only a small region beyond the QM system**
    * **Base your calculations on many MD simulation snapshots**
    * **Don't just minimise the energy, calculate the free energy (include dynamics)**



Choosing protonation states
---------------------------

How should one choose protonation states of aminoacids inside the
protein?  In the “methods” section of a publication one sometimes
finds: “Protonation states were chosen based on pKa values, except
Asp10, Glu43 and His35, which were protonated”, but without a decent
explanation.

* **pKa values inside the protein are far away from the solution, advanced methods needed to account for it. Most enzymologists use** `PROPKA <https://github.com/jensengroup/propka>`_ (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2876>`_) **.**


--------------------------------
Model preparation and validation
--------------------------------


Should one perform QM/MM calculation of a fully solvated protein in a periodic box of waters, droplet, implicit solvent or combination? If so, equilibrate at MM level, or not?


* **It is good to compare reactions energetics between gas phase, water and enzymes. One should see a clear catalytic effect in the protein environment.** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1092>`_)


* **When it comes to comparison with experimental data it is better to use explicit solvent.** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3257>`_)


* **QM cluster and QM/MM energies depend on the size of the QM system.** (`listen to Ulf Ryde on YouTube <https://youtu.be/aQdjC-W9Wy4?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2607>`_)


* **Recommendations and the "Big QM" cluster approach to get stable energies while increasing the number of QM atoms to 800-1000 atoms** (`listen to Ulf Ryde on YouTube <https://youtu.be/aQdjC-W9Wy4?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2905>`_) **:**
    * **Include in the QM regions neutral groups up to 4-5 Å away from the minimal QM system consisting of the active site**
    * **Include ALL the charged groups that are not on the surface of the protein (i.e. buried in the protein)**
    * **Move the “junction” atoms 2 residues +caps aways from the active site / minimal QM system**

      
* **QM/MM structures (and energies) are much more stable than QM-cluster structures (and energies) while increasing the QM size: smaller QM parts can be used with QM/MM models** (`listen to Ulf Ryde on YouTube <https://youtu.be/aQdjC-W9Wy4?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3061>`_)

|

====================================================
QM/MM modelling / simulation protocol and validation
====================================================

-------------------------------------
How best to choose a level of theory?
-------------------------------------

* **Projector-based embedding schemes are beneficial in getting consistent results for different DFT functionals** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1838>`_)


* **Coupled Cluster (CC) methods are accurate but slow, MP2-based methods practically are more reliable** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3300>`_)


* **First always refer to relevant literature for your specific problem** (`listen to Maria Khrenova on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=681>`_)


* **If nobody has studied your biological system yet, focus on the specific chemistry involved in the phenomenon you want to study and start looking at the levels of theory employed to study it described in the literature** (`listen to Maria Khrenova on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2151>`_)


* **Which Hamiltonian to choose** (`listen to Maria João Ramos on YouTube <https://youtu.be/XIHMcR_tR7E?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=639>`_) :
    * **From the literature you have to infer the approach that describes all the energy contributions involved in the phenomenon you are investigating**
    * **Consider the availability of software that implements it**
    * **Consider the availability of computational resources necessary to run it**


* **Recommendations for the method to choose** (`listen to Ulf Ryde on YouTube <https://youtu.be/aQdjC-W9Wy4?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3947>`_) **:**
    * **If structure properties are the aim: pure DFT + dispersion corrections + small basis sets**
    * **If energies are the aim: single points with larger basis sets**
    * **Test pure and hybrid functionals: if the results are not comparable then calibrate your energies with higher level of theory (e.g. CCSD(T) for closed shell case or DMRG-PT2 for an open shell system)**

      
* **The usage of polarization functions in the basis set used to describe the QM region is essential. Use at least a DZP basis set.** (`listen to Maria Khrenova on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2903>`_)


* **Limitation in the use of diffuse functions in a QM/MM setup (if required they can be employed to describe atoms only in the middle of the QM box)** (`listen to Maria Khrenova on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2864>`_)


* **Mixing the levels of theory to draw a conclusion could be dangerous** (`listen to Maria Khrenova on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1752>`_)
  
---------------------------------------------
How best to choose a suitable DFT functional?
---------------------------------------------

* **Benchmark your DFT functionals before embarking on expensive QM/MM calculations, modelling at least the relevant part of your system (e.g. active site) and taking a high-level theoretical method (e.g. CCSD(T)/CBS) as a reference for the energies.** (`listen to Maria João Ramos on YouTube <https://youtu.be/XIHMcR_tR7E?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=832>`_)

* **Do not use QM/MM (or any other method or computational tool) as a black box: if there are disagreements you need to explain them. Different functionals, inclusion of dispersion correction etc. can yield not only quantitavely but also qualitatively different results, and this can be to how well the approach captures - or fails to capture - key underlying chemistry.** (`listen to Maria Khrenova present 5 examples on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1181>`_)
  
* **Test the functional against experimental results** (`listen to Maria Khrenova on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3065>`_)

  
* **DFT is not a systematically improvable method, sometimes you need to go beyond that, DFT often gives too low barriers** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=843>`_)


* **Dispersion corrections for DFT can often improve results** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2527>`_)


* **Check if your DFT functional preserves reactants configurations well: conformation of the substrate, H-bonds, etc.** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1041>`_)


---------------------------------------------
How best to choose a suitable QM region size?
---------------------------------------------

* **QM/MM works very well for enzymes even with relatively small QM region sizes.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1730>`_)


* **One should focus more on the quality of QM treatment, rather than QM region size. QM/MM often doesn't converge with respect to the size of the QM part.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2093>`_)


* **One could check the effect of individual residues on catalysis by calculating their individual contributions into transition state stabilization.** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2525>`_)


* **To choose the initial QM region size it is good to check interactions at the MM level.** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2624>`_)


* **Number of QM atoms typically depends on (limited by) the available computational resources, in practice ~150 QM atoms should be reachable.** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3434>`_)


* **One protocol to find a suitable QM region, aiming to include all important effects in the QM region** (`listen to Ulf Ryde on YouTube <https://youtu.be/aQdjC-W9Wy4?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3897>`_) **:**
    * **Start with a rather small QM region and perform a QM/MM optimization with fixed surrounding**
    * **Repeat it with free surroundings**
    * **If there is a large difference between the results in the two previous points, then increase the QM size and repeat the cycle from step one**


-----------------------------------------------------------------------
How best to choose the valence saturation scheme at the QM-MM boundary?
-----------------------------------------------------------------------

* **Typically Link-atoms should be used.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=401>`_)


* **Put the QM-MM boundary preferably on an sp3 hybridised C atom (cut C-C bonds).** (`listen to Maria Khrenova on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2626>`_)


--------------------------------------------------------------------------------
Is the electrostatic coupling always the best compromise for the QM/MM coupling?
--------------------------------------------------------------------------------

* **Importance of the electrostatic embedding scheme. The mechanical scheme often brings underestimations and wrong results.** (`listen to Maria Khrenova on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2297>`_)


* **QM-MM interactions should also involve LJ-VdW interactions, not only electrostatics.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=348>`_)


------------------------------------------------------------------
What kind of systematic basis set benchmarking should one perform?
------------------------------------------------------------------

* **Check for convergence with respect to basis set size.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3986>`_)


---------------------------------------------------------------------------------------------
How to choose a good collective variable (reaction coordinate) and perform effective sampling
---------------------------------------------------------------------------------------------

How to choose a good collective variable (reaction coordinate) and perform effective sampling: chemical intuition versus unbiased and automated approaches.

* **It is good to start from the reaction in solution and/or cluster model within polarizable continuum to test possible reaction pathways.** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=710>`_)


* **Collective variable should include all bonds that form or break during the reaction:** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2272>`_)


* **Low-energy vibrations and variables in most cases do not affect free energy landscapes dramatically, so they could be excluded from collective variables.** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3529>`_)


* **Position restraints may affect your final results significantly and should not be used in actual profile simulations.** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3775>`_)


-------------------------------------------
Long-range electrostatics: truncate or not?
-------------------------------------------

* **Evaluate how important long-range (7-20 Å) interactions are for the problem you are dealing with and choose the right model to describe them.** (`listen to Maria João Ramos on YouTube <https://youtu.be/XIHMcR_tR7E?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1348>`_)

|
  
=============================================================
Validation, analysis and interpretation of simulation results
=============================================================

-------------------------------------------------------------------------
High-level QM with limited or no sampling, or low-level QM with sampling?
-------------------------------------------------------------------------

* **On a high level you could get only enthalpies.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1160>`_)


* **Do analysis of the free energy landscape: extract and analyse all reactants, products, transition state geometries and characterise them.** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2446>`_)


-------------------------------------------------------------------------------
How can we combine low-level QM for sampling with high-level QM for energetics?
-------------------------------------------------------------------------------

* **Static methods could give a good estimate of enthalpies. To get activation free energies one needs to consider dynamics. You could estimate entropic factors with lower-level QM/MM dynamics.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1075>`_)

* **To account for temperature effects one should consider QM/MD based free energy methods, like metadynamics.** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1794>`_)

---------------------------------------------------------------
How to check validity and convergence of the reaction pathways?
---------------------------------------------------------------

* **The reaction space (in particular for enzymatic reactions) can be very complex: the more complex the Hamiltonian, the more difficult may be the reaction mechanism (i.e. more steps involved).** (`listen to Maria João Ramos on YouTube <https://youtu.be/XIHMcR_tR7E?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1717>`_)

  
* **Sampling time is important, be careful and check free energy profile convergence.** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3331>`_)

  
* **In enzymes, the transition state should be stabilized by the environment. Check for stabilization factors relative to the solvent.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1286>`_)

  
* **Reactive conformation in the enzyme has the lowest activation energy, however they are often not the most populated.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2721>`_)


* **Sometimes Michaelis complex in the enzymes can have a distorted geometry, which is not energetically most favorable but resembles the transition state, thus lowering activation energy.** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=606>`_)


* **Beware of sampling in QM/MM, convergence could take a long time, especially with respect to the orientation of the reactants.** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=962>`_)

  
* **Michaelis complex formation constant is important for very fast enzymes.** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=3360>`_)


* **Note that the criteria you identify to describe when a reaction can take place may depend on the QM approach you chose.** (`listen to Maria Khrenova on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=885>`_)


------------------------------------
How to validate the final result(s)?
------------------------------------

* **Try to make your prediction for well-studied systems and validate against experimental data.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=458>`_)


* **All intermediates in the reaction should be lower in energy than experimental activation energy.** (`listen to Adrian Mulholland on YouTube <https://youtu.be/8PGHcNKOLqY?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=729>`_)


* **Beware of improper models, some  effects do not contribute to the catalysis.** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=568>`_)


* **Tunneling of protons is perfectly fine with transition state theory if it is properly accounted for.** (`listen to Janez Mavri on YouTube <https://youtu.be/GjRTQ5Q13qg?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1354>`_)


* **If you want to design an enzyme inhibitor, take into account that it should resemble both reactants state of the substrate as well as adapt for a binding site charge distribution well.** (`listen to Carme Rovira on YouTube <https://youtu.be/mojq6K6N7UM?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=2894>`_)


* **If you have to verify a hypothesis, first make direct tests, but then, make predictions with that hypothesis and verify them.** (`listen to Maria Khrenova on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1274>`_)


* **Test robustness of the results by benchmarking the level of theory (functionals + basis sets).** (`listen to Maria Khrenova on YouTube <https://youtu.be/uP1px6Yul2s?list=PLzLqYW5ci-2d-wolQ9CpE4akorB3naRso&t=1063>`_)

  
