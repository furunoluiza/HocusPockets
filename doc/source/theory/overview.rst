Theoretical Overview
====================

HocusPockets integrates three complementary approaches to identify and
rank druggable cavities in protein structures.

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

3. Cavity scoring
-----------------

Each detected cavity is characterized by:

* **Geometric descriptors**: volume (A³), surface area (A²), diameter,
  radius, minimal thickness, enclosing cylinder radius
* **Physicochemical descriptors**: average hydrophobicity (Wimley-White
  scale), aromatic residue frequency, tyrosine OH frequency
* **Reactivity descriptors**: mean Fukui and Fukushima values across
  cavity atoms

These descriptors are combined into a composite **Score** used for
ranking potential binding sites.
