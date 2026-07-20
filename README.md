<p align="center"> <img src="doc/source/_static/HocusPockets_logo.png" alt="HocusPockets" /> </p>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![Documentation](https://img.shields.io/badge/docs-Sphinx-purple.svg)](https://furunoluiza.github.io/HocusPockets/)
[![DOI](https://img.shields.io/badge/doi-10.1101/2024.01.01.575590-blue)]

HocusPockets (**Hol**istic **C**omputational **U**tility for **S**coring **Pockets**) is a tool designed to prioritize biological targets within the proteome of microorganisms based on the druggability of interaction sites and their chemical reactivity, utilizing descriptors such as Fukui indices to evaluate potential binding sites.

## Prerequisites

Before using HocusPockets, ensure the following packages and external programs are installed:

### Python packages

> * pandas  
> * numpy  
> * scipy
> * matplotlib


### External software

> * [Fpocket](https://github.com/Discngine/fpocket/tree/master)  
> * [ConCavity](https://compbio.cs.princeton.edu/concavity/)  
> * [PRIMoRDiA](https://github.com/igorChem/PRIMoRDiA1.0v)  
> * [MOPAC](https://openmopac.net/)
    
## Usage

You can use HocusPockets with either a .pdb or .cif file format:
```bash
hocuspockets -i file.pdb -o output_path
```

For a full list of options, run:
```bash
hocuspockets -h
```

For more information, please refer to the [documentation](https://furunoluiza.github.io/HocusPockets/).