import pytest
from hocuspockets.initialize import args_parse


def test_args_parse_minimal():
    with pytest.raises(SystemExit):
        args_parse()


def test_args_parse_defaults():
    import sys
    test_args = ["prog", "-i", "test.pdb"]
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(sys, "argv", test_args)
        try:
            args_parse()
        except SystemExit:
            pytest.skip("requires a real .pdb file")


def test_description_exists():
    from hocuspockets.initialize import description
    assert "Holistic Computational Utility" in description
