# Description
###############################################################################
"""
This module provides utility functions for parsing and manipulating PDB
(Protein Data Bank) files used throughout the HocusPockets pipeline.

It centralizes common PDB parsing operations such as extracting residue IDs,
atom coordinates, residue centers, and parsing B-factor values. Using these
utilities reduces code duplication across the cavity analysis modules.
"""

# Imports
###############################################################################
import logging
from typing import Dict, List, Optional, Set, Tuple
import numpy as np
import numpy.typing as npt

logger = logging.getLogger(__name__)

# License
###############################################################################
"""
    HocusPockets: Holistic Computacional Utility for Scoring Pockets

    Authors: Machado, L.F.; Torres, P.H.M.; Spelta, G.I.; Oliveira Filho, J.A., Pascutti, P.G.

    [Federal University of Rio de Janeiro]

    Contact info:
    Carlos Chagas Filho Institute of Biophysics
    LMDM - Laboratory for Molecular Modeling and Dynamics
    Carlos Chagas Filho Av. 373 - Rio de Janeiro - RJ - Brazil
    Health Sciences Centre - block G1-19
    E-mail address: luizafuruno@ufrj.br
"""

PDB_COLUMNS: Dict[str, Tuple[int, int]] = {
    "record": (0, 6),
    "atom_serial": (6, 11),
    "atom_name": (12, 16),
    "residue_name": (17, 20),
    "chain": (21, 22),
    "residue_seq": (22, 26),
    "x": (30, 38),
    "y": (38, 46),
    "z": (46, 54),
    "occupancy": (54, 60),
    "bfactor": (60, 66),
}


def parse_pdb_lines(file_path: str, record_type: str = "ATOM") -> List[str]:
    """
    Reads a PDB file and returns lines matching the specified record type.

    Parameters
    ----------
        file_path : str
            Path to the PDB file.
        record_type : str, optional
            Record type to filter (e.g. "ATOM", "HETATM"). Default is "ATOM".

    Return
    ------
        list
            List of lines from the PDB file that start with the given record type.
    """
    lines: List[str] = []
    try:
        with open(file_path, "r") as f:
            for line in f:
                if line.startswith(record_type):
                    lines.append(line)
    except FileNotFoundError:
        logger.warning("PDB file not found: %s", file_path)
        return lines
    logger.debug("Parsed %d %s lines from %s", len(lines), record_type, file_path)
    return lines


def extract_coordinates(line: str) -> npt.NDArray[np.float64]:
    """
    Extracts the (x, y, z) coordinates from a PDB ATOM/HETATM line.

    Parameters
    ----------
        line : str
            A PDB line containing atomic coordinates in columns 30-54.

    Return
    ------
        numpy.ndarray
            1D array with three float values [x, y, z].
    """
    return np.array(
        [
            float(line[30:38].strip()),
            float(line[38:46].strip()),
            float(line[46:54].strip()),
        ]
    )


def extract_residue_id(line: str) -> str:
    """
    Extracts a unique residue identifier from a PDB line.

    Parameters
    ----------
        line : str
            A PDB ATOM line.

    Return
    ------
        str
            Residue ID in the format "chain:seqResName"
            (e.g. "A:123ALA").
    """
    chain_id = line[21].strip()
    res_seq = line[22:26].strip()
    res_name = line[17:20].strip()
    return f"{chain_id}:{res_seq}{res_name}"


def extract_atom_key(line: str) -> str:
    """
    Extracts a unique atom key from a PDB line.

    Parameters
    ----------
        line : str
            A PDB ATOM line.

    Return
    ------
        str
            Atom key in the format "chain:seqResName:atomName"
            (e.g. "A:123ALA:CA").
    """
    chain_id = line[21].strip()
    res_seq = line[22:26].strip()
    res_name = line[17:20].strip()
    atom_name = line[12:16].strip()
    return f"{chain_id}:{res_seq}{res_name}:{atom_name}"


def extract_residue_ids(file_path: str) -> Set[str]:
    """
    Extracts all unique residue identifiers from a PDB file.

    Parameters
    ----------
        file_path : str
            Path to the PDB file.

    Return
    ------
        set
            Set of residue IDs in "chain:seqResName" format.
    """
    residues: Set[str] = set()
    for line in parse_pdb_lines(file_path):
        residues.add(extract_residue_id(line))
    return residues


def calculate_residue_centers(lines: List[str]) -> Dict[str, npt.NDArray[np.float64]]:
    """
    Calculates the geometric center (centroid) of each residue from ATOM lines.

    Parameters
    ----------
        lines : list
            List of PDB ATOM lines.

    Return
    ------
        dict
            Residue ID (str) -> centroid coordinates (numpy.ndarray of shape (3,)).
    """
    residues: Dict[str, List[npt.NDArray[np.float64]]] = {}
    for line in lines:
        if line.startswith("ATOM"):
            residue_id = extract_residue_id(line)
            coords = extract_coordinates(line)
            if residue_id not in residues:
                residues[residue_id] = []
            residues[residue_id].append(coords)
    return {
        residue_id: np.mean(coords_list, axis=0)
        for residue_id, coords_list in residues.items()
    }


def extract_atom_coordinates(lines: List[str]) -> npt.NDArray[np.float64]:
    """
    Extracts all (x, y, z) coordinates from ATOM and HETATM lines.

    Parameters
    ----------
        lines : list
            List of PDB lines.

    Return
    ------
        numpy.ndarray
            Array of shape (N, 3) with atomic coordinates,
            or an empty array if no valid coordinates are found.
    """
    coords: List[List[float]] = []
    skipped = 0
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                c = extract_coordinates(line).tolist()
                coords.append(c)
            except (ValueError, IndexError):
                skipped += 1
                continue
    if skipped:
        logger.debug("Skipped %d lines with unparseable coordinates", skipped)
    if not coords:
        logger.debug("No valid coordinates extracted from %d lines", len(lines))
        return np.empty((0, 3))
    logger.debug("Extracted %d coordinates", len(coords))
    return np.array(coords)


def calculate_geometric_center(lines: List[str]) -> npt.NDArray[np.float64]:
    """
    Calculates the geometric center (centroid) from a set of PDB atomic lines.

    Parameters
    ----------
        lines : list
            List of PDB ATOM/HETATM lines.

    Return
    ------
        numpy.ndarray
            1D array with (x, y, z) coordinates of the centroid,
            or (0.0, 0.0, 0.0) if no coordinates are present.
    """
    coords = extract_atom_coordinates(lines)
    if len(coords) == 0:
        return np.array([0.0, 0.0, 0.0])
    return np.mean(coords, axis=0)


def build_residue_index_map(pdb_path: str) -> Dict[int, str]:
    """
    Builds a sequential index map for residues in a PDB file.

    Parameters
    ----------
        pdb_path : str
            Path to the PDB file.

    Return
    ------
        dict
            1-based sequential index (int) -> residue ID (str)
            in "chain:seqResName" format.
    """
    residue_map: Dict[int, str] = {}
    seen: Set[str] = set()
    index = 1
    for line in parse_pdb_lines(pdb_path):
        res_id = extract_residue_id(line)
        if res_id not in seen:
            seen.add(res_id)
            residue_map[index] = res_id
            index += 1
    return residue_map


def parse_bfactor_values(
    pdb_file: str, residue_map: Dict[int, str]
) -> Dict[str, float]:
    """
    Parses B-factor values from a PDB file, mapping atoms to their B-factors.

    Parameters
    ----------
        pdb_file : str
            Path to the PDB file.
        residue_map : dict
            Dictionary mapping sequential residue indices (int) to residue IDs (str),
            as returned by build_residue_index_map().

    Return
    ------
        dict
            Atom key (str) in "chain:seqResName:atomName" format -> B-factor value (float).
    """
    atom_values: Dict[str, float] = {}
    current_res_number: Optional[int] = None
    res_id: Optional[str] = None
    seq_index = 0
    skipped = 0
    mapped = 0

    try:
        fh = open(pdb_file, "r")
    except FileNotFoundError:
        logger.warning("B-factor PDB file not found: %s", pdb_file)
        return atom_values

    with fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue

            try:
                res_number = int(line[22:26].strip())
                atom_name = line[12:16].strip()
                bfactor = float(line[60:66].strip())
            except ValueError:
                skipped += 1
                continue

            if res_number != current_res_number:
                seq_index += 1
                res_id = residue_map.get(seq_index)
                current_res_number = res_number

            if not res_id:
                skipped += 1
                continue

            atom_key = f"{res_id}:{atom_name}"
            atom_values[atom_key] = bfactor
            mapped += 1

    logger.debug(
        "Parsed %d B-factor values from %s (%d lines skipped)",
        mapped,
        pdb_file,
        skipped,
    )
    return atom_values
