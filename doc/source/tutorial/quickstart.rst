Quickstart Tutorial

This tutorial walks you through a complete HocusPockets analysis using a
protein structure in PDB format.

Step 1: Prepare your input
--------------------------

Obtain a protein structure file in PDB or mmCIF format from the
`Protein Data Bank <https://www.rcsb.org/>`_. For example:

.. code-block:: bash

    wget https://files.rcsb.org/download/1a30.pdb

PDB files need to be preprocessed to remove water molecules, ligands and protonated residues. 
You can use the pdb2pqr tool for this purpose.

Step 2: Run HocusPockets
-------------------------

.. code-block:: bash

    hocuspockets -i 1a30.pdb -o results/

Step 3: Inspect the output
---------------------------

After execution, the output directory will contain:

``cavities/``
    Subdirectories for each detection method (fpocket, ligsite, surfnet,
    pocketfinder, fukui) with grid and residue PDB files.

``reactivity/``
    MOPAC and PRIMoRDiA outputs, including Fukui and Fukushima descriptor
    PDB files.

``<protein>_cavity_scores.csv``
    A ranked table of all detected cavities with geometric and
    physicochemical descriptors.

``<protein>_cavity_overlap.csv``
    Pairwise residue overlap between cavities from different methods.

Step 4: Interpret scores
-------------------------

Cavities are ranked by a composite ``Score`` that blends:

* Geometric properties (volume, surface, diameter)
* Hydrophobicity and aromatic residue frequency
* Quantum reactivity (Fukui and Fukushima indices)

Higher scores indicate more promising druggable cavities.
Quickstart Tutorial

This tutorial walks you through a complete HocusPockets analysis using a
protein structure in PDB format.

Step 1: Prepare your input
--------------------------

Obtain a protein structure file in PDB or mmCIF format from the
`Protein Data Bank <https://www.rcsb.org/>`_. For example:

.. code-block:: bash

    wget https://files.rcsb.org/download/1a30.pdb

Step 2: Run HocusPockets
-------------------------

.. code-block:: bash

    hocuspockets -i 1a30.pdb -o results/

Step 3: Inspect the output
---------------------------

After execution, the output directory will contain:

``cavities/``
    Subdirectories for each detection method (fpocket, ligsite, surfnet,
    pocketfinder) with grid and residue PDB files.

``cavities/ranked_cavities/``
    The five most promising cavities, ranked according to their
    Fukushima sum.

``reactivity/``
    MOPAC and PRIMoRDiA outputs, including Fukui and Fukushima descriptor
    PDB files.

``<protein>_cavity_scores.csv``
    A ranked table of all detected cavities with geometric and
    physicochemical descriptors.

``<protein>_cavity_overlap.csv``
    Pairwise residue overlap between cavities from different methods.

Step 4: Cavity characterization
-------------------------

HocusPockets characterizes each detected cavity using a combination of
geometric, physicochemical, and electronic-reactivity descriptors.

* Geometric descriptors: volume, surface area, diameter, radius,
    smallest dimension, and cylinder radius.
* Physicochemical descriptors: hydrophobicity and aromatic residue
    frequency, including the frequency of tyrosine residues.
* Electronic-reactivity descriptors: Fukui and Fukushima descriptors
    calculated for residues associated with each cavity.

These descriptors provide a multidimensional characterization of the
structural and chemical properties of each predicted binding cavity.

Step 5: Cavity prioritization
-------------------------
Cavities are prioritized according to their **Fukushima Sum**, which
represents the cumulative Fukushima reactivity descriptor associated
with the residues comprising each cavity.

The five cavities with the highest Fukushima sums are selected as the
most promising cavities and stored in:

``cavities/ranked_cavities/``

A higher Fukushima sum indicates greater cumulative electronic
reactivity within the cavity and is used as the primary criterion for
cavity prioritization.
