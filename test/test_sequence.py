"""
module to test sequence.py
"""
import pytest
from seq_tools.sequence import to_dna, to_rna, get_reverse_complement

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
