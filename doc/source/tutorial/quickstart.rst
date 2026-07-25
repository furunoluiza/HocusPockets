Quickstart Tutorial
===================

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
