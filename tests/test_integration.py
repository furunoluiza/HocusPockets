import os
import pytest
from hocuspockets import pdb_utils

PDBBIND_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "database_pdbbind")


@pytest.mark.skipif(not os.path.isdir(PDBBIND_DIR), reason="PDBbind database not available")
def test_real_pdb_parsing():
    pdb_path = os.path.join(PDBBIND_DIR, "1a30.pdb")
    assert os.path.isfile(pdb_path)
    residues = pdb_utils.extract_residue_ids(pdb_path)
    assert len(residues) > 0
    coords = pdb_utils.extract_atom_coordinates(pdb_utils.parse_pdb_lines(pdb_path))
    assert coords.shape[0] > 0
    assert coords.shape[1] == 3


@pytest.mark.skipif(not os.path.isdir(PDBBIND_DIR), reason="PDBbind database not available")
def test_real_pdb_residue_centers():
    pdb_path = os.path.join(PDBBIND_DIR, "1a30.pdb")
    centers = pdb_utils.calculate_residue_centers(pdb_utils.parse_pdb_lines(pdb_path))
    assert len(centers) > 0
