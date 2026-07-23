# HocusPockets — Agent guide

## Quick start

```bash
pip install -e ".[dev]"                 # editable install + dev deps (ruff, pytest)
hocuspockets -i file.pdb -o results/    # run the pipeline
hocuspockets-install                    # download/install external tools via conda
```

## Non-Python requirements (must be on PATH)

The pipeline exits immediately if any is missing:

- `fpocket` — `conda install -c conda-forge fpocket`
- `concavity` — download + `make` from https://compbio.cs.princeton.edu/concavity/
- `PRIMoRDiA_1.25v_LINUX64` — binary from https://github.com/igorChem/PRIMoRDiA1.0v
- `mopac` — `conda install -c conda-forge mopac`

`hocuspockets-install` automates these via a `conda run -n HocusPockets` environment.

## CLI

```
hocuspockets -i PATH [-o PATH] [--distance FLOAT] [--similarity FLOAT]
              [--fpocket-score FLOAT] [--eband-limit INT]
              [--parallel] [--n-jobs INT] [--verbose]
```

- `-i` accepts a single `.pdb` file or a directory of `.pdb` files **only**. The README mentions `.cif` but the code rejects anything not `.pdb`.
- **Directory input gotcha**: if the directory contains *any* non-`.pdb` file, the pipeline exits with an error (`initialize.py:259-261`). Only pure `.pdb` directories are accepted.
- `--similarity` must be in `[0.01, 1.0]` (step 0.01; enforced via argparse `choices`).
- `--verbose` writes `verbose.log` to the output directory.
- `--parallel` uses `ThreadPoolExecutor` (not multiprocessing).

## Pipeline (execution order)

1. Tool availability checks → `fpocket`, `concavity`, `PRIMoRDiA`, `mopac`
2. `Reactivity.find_reactivity()` → MOPAC (PM7, 1SCF, MOZYME) → PRIMoRDiA → Fukui/Fukushima descriptors
3. `CavityFinder.find_cavities()`: fpocket → ConCavity (`pocketfinder`, `ligsite`, `surfnet`) → catalytic cavity via top-5% Fukui → cavity overlap Jaccard comparison → physicochemical scoring CSV

## Project structure

```
hocuspockets/              # single package
  __main__.py              # entry point → main()
  initialize.py            # CLI args, tool checks, color dict
  reactivity.py            # Reactivity class (MOPAC + PRIMoRDiA)
  cavity_finder.py         # CavityFinder, FpocketFinder, ConCavityFinder
  scoring_pockets.py       # ScoringPockets (physicochemical descriptors, CSV output)
  catalytic_cav.py         # CatalyticCav (Fukui/Fukushima cavity averages)
  install_external.py      # download/install helper
analyse_results/           # standalone analysis scripts (not part of package)
database_pdbbind/          # 270 PDB files (potential test fixtures)
doc/                       # Sphinx docs: cd doc && make html
```

Pipeline output hierarchy: `{output_path}/{protein_name}/cavities/{fpocket,ligsite,surfnet,pocketfinder,fukui}/` plus `{protein_name}_cavity_scores.csv` and `{protein_name}_cavity_overlap.csv`.

## Development

- **Lint**: `ruff check hocuspockets/` (dev dep)
- **Tests**: none exist. Existing PDB fixtures in `database_pdbbind/` could be used.
- **Docs**: `cd doc && make html` (Sphinx + `sphinx-rtd-theme`, requires `pip install -e ".[docs]"`)
- **Build**: setuptools via `pyproject.toml`; entry points `hocuspockets` and `hocuspockets-install`.
- **MOPAC caveat**: uses PM7 + 1SCF + MOZYME + COSMO (eps=78.4). May fail on unusual residues or fail to converge.
