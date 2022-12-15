"""
module to test sequence.py
"""

from seq_tools.sequence import (
    to_dna,
    to_rna,
    get_reverse_complement,
    get_molecular_weight,
    get_max_stretch,
)


def test_to_dna():
    """
    test to_dna function
    """
    assert to_dna("AUG") == "ATG"
    assert to_dna("AUGAUGAUG") == "ATGATGATG"
    assert to_dna("AUGAUGAUGAUG") == "ATGATGATGATG"


def test_to_rna():
    """
    test to_rna function
    """
    assert to_rna("ATG") == "AUG"
    assert to_rna("ATGATGATG") == "AUGAUGAUG"
    assert to_rna("ATGATGATGATG") == "AUGAUGAUGAUG"


def test_get_reverse_complement():
    """
    test get_reverse_complement function
    """
    assert get_reverse_complement("AUG", "RNA") == "CAU"
    assert get_reverse_complement("CAU", "RNA") == "AUG"
    assert get_reverse_complement("ATGATGATG") == "CATCATCAT"
    assert get_reverse_complement("ATGATGATGATG") == "CATCATCATCAT"


def test_get_molecular_weight():
    """
    test get_molecular_weight function
    """
    assert get_molecular_weight("AUG", "RNA") == 1034.6
    assert get_molecular_weight("CAU", "RNA") == 994.5999999999999
    assert get_molecular_weight("ATGATGATG") == 3001.7999999999997
    assert (
        get_molecular_weight("ATGATGATGATG", double_stranded=True) == 7844.799999999998
    )


def test_get_max_stretch():
    """
    test get_max_stretch function
    """
    assert get_max_stretch("AUG") == 1
    assert get_max_stretch("CAU") == 1
    assert get_max_stretch("AGGA") == 2
    assert get_max_stretch("AAAACCC") == 4
