<p align="center"> <img src="doc/source/_static/HocusPockets_logo.png" height="110" alt="HocusPockets" /> </p>
HocusPockets is a tool designed to prioritize biological targets within the proteome of microorganisms based on the druggability of interaction sites and their chemical reactivity, utilizing descriptors such as Fukui indices to evaluate potential binding sites.

## Prerequisites

Before using HocusPockets, ensure the following packages and external programs are installed:

### Python packages

> * pandas  
> * numpy  
> * scipy
> * matplotlib


### External software

> * [fpocket](https://github.com/Discngine/fpocket/tree/master)  
> * [ConCavity](https://compbio.cs.princeton.edu/concavity/)  
> * [PRIMoRDiA](https://github.com/igorChem/PRIMoRDiA1.0v)  
    
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