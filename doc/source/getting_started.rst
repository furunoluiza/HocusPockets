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

If you want to install external software automatically, you can use the
following command:

.. code-block:: bash

    hocuspockets-install

The automatic installer builds PRIMoRDiA 1.50 with CMake and installs it under
your home directory. PRIMoRDiA 1.50 requires CMake, a C++17 compiler with
OpenMP, Eigen, Boost.Iostreams, zlib, and Make. On Ubuntu/Debian, the native
build dependencies can be installed with:

.. code-block:: bash

    sudo apt update
    sudo apt install build-essential cmake libeigen3-dev libboost-iostreams-dev zlib1g-dev

Open a new terminal and activate the conda environment, or source the
changes in the current terminal:

.. code-block:: bash

    source ~/.bashrc

Quickstart
----------

Display usage options:

.. code-block:: bash

    hocuspockets --help

You can use HocusPockets with either a .pdb or .cif file format:

.. code-block:: bash

    hocuspockets -i file.pdb -o output_path
