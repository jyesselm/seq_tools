import pytest
from seq_tools.sequence import *

def test_get_nucleic_acid_type():
    seq = "ACGT"
    assert get_nucleic_acid_type(seq) == "DNA"
    seq = "ACGU"
    assert get_nucleic_acid_type(seq) == "RNA"
    seq = "AAA"
    assert get_nucleic_acid_type(seq) == "DNA"


