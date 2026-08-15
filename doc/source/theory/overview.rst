Theoretical Overview

HocusPockets integrates three complementary approaches to identify and rank druggable cavities in protein structures.

1. Geometric cavity detection
-----------------------------

**Fpocket** uses alpha-spheres to detect pockets on the protein surface
based purely on geometry (van der Waal radii and solvent accessibility).

**ConCavity** extends pocket detection by incorporating evolutionary
conservation. It supports three grid methods:

* ``pocketfinder`` — knowledge-based potentials
* ``ligsite`` — geometric detection of clefts
* ``surfnet`` — solvent-accessible surface mapping

2. Quantum reactivity descriptors
----------------------------------

**MOPAC** performs semi-empirical quantum chemistry (PM7 Hamiltonian) to
compute molecular orbitals. **PRIMoRDiA** then calculates:

**Fukui indices (f⁺/f⁻)**
    Measure of electrophilic / nucleophilic reactivity at each atomic
    site. High values indicate atoms likely to participate in chemical
    reactions.

**Fukushima indices**
    Derived from orbital-weighted Fukui potentials, they capture the
    regioselective reactivity relevant for covalent drug design.
Theoretical Overview

HocusPockets integrates complementary structure-based and quantum-chemical
approaches to detect, characterize, and prioritize potential protein
binding cavities.

1. Geometric cavity detection
-----------------------------

**Fpocket** uses alpha-spheres to detect pockets on the protein surface
based purely on geometry (van der Waal radii and solvent accessibility).

**ConCavity** extends pocket detection by incorporating evolutionary
conservation. It supports three grid methods:

* ``pocketfinder`` — knowledge-based potentials
* ``ligsite`` — geometric detection of clefts
* ``surfnet`` — solvent-accessible surface mapping

2. Quantum reactivity descriptors
----------------------------------

HocusPockets incorporates quantum-chemical calculations to characterize
the electronic reactivity of residues associated with detected cavities.

**MOPAC** performs semi-empirical quantum-chemical calculations using the
PM7 Hamiltonian to obtain the electronic structure of the molecular
system. **PRIMoRDiA** subsequently uses these results to calculate
electronic-reactivity descriptors.

**Fukui indices (f⁺/f⁻)**
    Measure of electrophilic / nucleophilic reactivity at each atomic
    site. High values indicate atoms likely to participate in chemical
    reactions.

**Fukushima indices**
    Derived from orbital-weighted Fukui potentials, they capture the
    regioselective reactivity relevant for covalent drug design.

3. Cavity characterization
-----------------

Each detected cavity is characterized using a combination of geometric,
physicochemical, and electronic-reactivity descriptors.

* **Geometric descriptors**: volume (A³), surface area (A²), diameter,
  radius, minimal thickness, and enclosing cylinder radius
* **Physicochemical descriptors**: average hydrophobicity (Wimley-White
  scale), and aromatic residue frequency
* **Reactivity descriptors**: mean Fukui, mean Fukushima, and cumulative Fukushima

These descriptors are combined into a composite **Score** used for
ranking potential binding sites.

4. Cavity prioritization
-----------------

Rather than combining the cavity descriptors into a composite score,
HocusPockets prioritizes cavities according to their cumulative
Fukushima reactivity.

For each cavity, the Fukushima sum represents the cumulative
Fukushima descriptor associated with its constituent residues. Cavities
with higher Fukushima sums are prioritized as having greater cumulative
electronic reactivity.

The five highest-ranked cavities are selected as the most promising
candidates and stored in the cavities/ranked_cavities/ directory.
