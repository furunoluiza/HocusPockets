from typing import Dict, Tuple

HYDROPHOBICITY_SCALE: Dict[str, float] = {
    "ALA": -0.17,
    "ARG": -0.81,
    "ASN": -0.42,
    "ASP": -1.23,
    "CYS": 0.24,
    "GLN": -2.02,
    "GLU": -0.58,
    "GLY": -0.01,
    "HIS": -0.96,
    "ILE": 0.31,
    "LEU": 0.56,
    "LYS": -0.99,
    "MET": 0.23,
    "PHE": 1.13,
    "PRO": -0.45,
    "SER": -0.13,
    "THR": -0.14,
    "TRP": 1.85,
    "TYR": 0.94,
    "VAL": -0.07,
    "MME": 0.23,
    "MSE": 0.23,
    "UNK": -0.17,
}

MOPAC_BINARY: str = "/opt/mopac/MOPAC2016.exe"
PRIMORDIA_BINARY: str = "PRIMoRDiA_1.25v_LINUX64"
GRID_METHODS: Tuple[str, ...] = ("pocketfinder", "ligsite", "surfnet")
CAVITY_PROGRAMS: Tuple[str, ...] = ("fpocket", "ligsite", "surfnet", "pocketfinder")
AROMATIC_RESIDUES: set = {"PHE", "TYR", "HIS", "TRP"}
DEFAULT_MOPAC_FLAGS: str = "PM7 1SCF ALLVECS VECTOR AUX LARGE MOZYME eps=78.4"
