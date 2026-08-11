import tempfile
import os
import pytest
from hocuspockets import pdb_utils
from conftest import SAMPLE_PDB_LINES


@pytest.fixture
def pdb_file():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.writelines(SAMPLE_PDB_LINES)
        tmp_path = f.name
    yield tmp_path
    os.unlink(tmp_path)


def test_extract_residue_ids(pdb_file):
    residues = pdb_utils.extract_residue_ids(pdb_file)
    assert len(residues) == 2
    assert "A:1ALA" in residues
    assert "B:2LYS" in residues


def test_extract_coordinates():
    line = SAMPLE_PDB_LINES[1]
    coords = pdb_utils.extract_coordinates(line)
    assert coords.shape == (3,)
    assert coords[0] == 7.500
    assert coords[1] == 6.891
    assert coords[2] == 6.379


def test_extract_residue_id():
    line = SAMPLE_PDB_LINES[1]
    res_id = pdb_utils.extract_residue_id(line)
    assert res_id == "A:1ALA"


def test_calculate_residue_centers():
    centers = pdb_utils.calculate_residue_centers(SAMPLE_PDB_LINES)
    assert len(centers) == 2
    assert "A:1ALA" in centers
    assert "B:2LYS" in centers
    ala_center = centers["A:1ALA"]
    assert len(ala_center) == 3


def test_extract_atom_coordinates():
    coords = pdb_utils.extract_atom_coordinates(SAMPLE_PDB_LINES)
    assert coords.shape[0] == 8
    assert coords.shape[1] == 3


def test_calculate_geometric_center():
    center = pdb_utils.calculate_geometric_center(SAMPLE_PDB_LINES)
    assert center.shape == (3,)


def test_build_residue_index_map(pdb_file):
    mapping = pdb_utils.build_residue_index_map(pdb_file)
    assert len(mapping) == 2
    assert mapping[1] == "A:1ALA"
    assert mapping[2] == "B:2LYS"


def test_parse_pdb_lines(pdb_file):
    lines = pdb_utils.parse_pdb_lines(pdb_file)
    assert len(lines) == 8
    lines_het = pdb_utils.parse_pdb_lines(pdb_file, record_type="HETATM")
    assert len(lines_het) == 0


def test_empty_residue_centers():
    centers = pdb_utils.calculate_residue_centers([])
    assert centers == {}


def test_empty_atom_coordinates():
    coords = pdb_utils.extract_atom_coordinates(["TER\n", "END\n"])
    assert coords.shape == (0, 3)
