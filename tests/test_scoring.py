import pytest
from hocuspockets.scoring_pockets import ScoringPockets


class DummyLogger:
    def info(self, msg):
        pass
    def warning(self, msg):
        pass
    def error(self, msg):
        pass


@pytest.fixture
def scorer():
    return ScoringPockets("", "", DummyLogger())


def test_score_cavity_basic(scorer):
    score = scorer.score_cavity(
        volume=100.0, surface=200.0, diameter=10.0, smallest_size=5.0,
        hydrophobicity_score=0.5, aromatic_freq=0.2, otyr_freq=0.1,
        radius=3.0, cylinder_radius=2.0,
    )
    assert isinstance(score, float)
    assert score > 0


def test_score_cavity_zero_values(scorer):
    score = scorer.score_cavity(
        volume=0.0, surface=0.0, diameter=0.0, smallest_size=0.0,
        hydrophobicity_score=0.0, aromatic_freq=0.0, otyr_freq=0.0,
        radius=0.0, cylinder_radius=0.0,
    )
    assert score == 0.0


def test_calculate_geometric_descriptors_insufficient_points(scorer):
    lines = [
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  ",
    ]
    result = scorer.calculate_geometric_descriptors(lines, "test.pdb")
    assert result == (0, 0, 0, 0, 0, 0)


def test_calculate_geometric_center_empty(scorer):
    center = scorer.calculate_geometric_center([])
    assert center == (0.0, 0.0, 0.0)
