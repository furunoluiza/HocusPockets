Getting Started
===============

Installation
------------

Clone the repository:

.. code-block:: bash

    git clone https://github.com/furunoluiza/HocusPockets

Enter the project folder:

.. code-block:: bash

    cd HocusPockets

Create a conda environment:

.. code-block:: bash

    conda create -n HocusPockets python=3.10
    conda activate HocusPockets

Install the package:

.. code-block:: bash

    pip install .

----

Quickstart
----------

Display usage options:

.. code-block:: bash

    HocusPockets --help

You can use HocusPockets with either a .pdb or .cif file format:
.. code-block:: bash

    HocusPockets -i file.pdb -o output_path

