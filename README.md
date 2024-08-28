# BINDiS: Best Interaction Network Druggability Score

BINDiS is a tool designed to prioritize biological targets within the proteome of microorganisms based on the druggability of interaction sites.

## Prerequisites

Before using BINDiS, ensure the following packages and external programs are installed:

### Python packages

> * os
> * sys
> * argparse
> * subprocess
> * shutil
> * textwrap

### External software
> * [fpocket](https://github.com/Discngine/fpocket/tree/master)

## Usage

You can use BINDiS with either a .pdb or .cif file format:
```bash
python3 ./main.py -f file.pdb -o output_path
```

For a full list of options, run:
```bash
python3 ./main.py -h
```
