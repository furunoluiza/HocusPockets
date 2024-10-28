# Hocus Pockus: Holistic Computacional Utility for Scoring Pockets

Hocus Pockus is a tool designed to prioritize biological targets within the proteome of microorganisms based on the druggability of interaction sites and their chemical reactivity, utilizing descriptors such as Fukui indices to evaluate potential binding sites.

## Prerequisites

Before using WhichCav, ensure the following packages and external programs are installed:

### Python packages

> * os
> * sys
> * argparse
> * subprocess
> * shutil

### External software
> * [fpocket](https://github.com/Discngine/fpocket/tree/master)

## Usage

You can use WhichCav with either a .pdb or .cif file format:
```bash
python3 ./WhichCav/ -f file.pdb -o output_path
```

For a full list of options, run:
```bash
python3 ./WhichCav/ -h
```
