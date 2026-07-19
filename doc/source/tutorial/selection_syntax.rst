Selection Syntax
================

HocusPockets accepts input in two formats:

* ``.pdb`` — Protein Data Bank format
* ``.cif`` / ``.mmcif`` — Macromolecular Crystallographic Information File

Input can be a single file or a directory containing multiple protein
structures:

.. code-block:: bash

    # Single file
    hocuspockets -i protein.pdb -o output/

    # Directory of PDB files
    hocuspockets -i ./pdbs/ -o output/

Options reference
-----------------

.. code-block:: text

    -i, --input PATH        Input file (.pdb/.cif) or directory (required)
    -o, --output PATH       Output directory (default: current dir)
    --distance VALUE        Max distance (A) between grid and residues (6.0)
    --similarity VALUE      Min residue similarity between pockets (0.5)
    --fpocket-score VALUE   Min druggability score for fpocket (0.5)
    --eband-limit VALUE     HOMO-LUMO gap limit in eV for orbital selection (3)
    --parallel              Enable parallel reactivity calculations
    --n-jobs VALUE          Number of CPU cores for parallel execution
    --verbose               Enable verbose logging to file
