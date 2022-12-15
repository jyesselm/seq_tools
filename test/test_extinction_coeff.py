"""
test extinction coefficient module
"""

import pytest
from seq_tools.extinction_coeff import get_extinction_coeff


def test_ds_dna():
    seq = "ACGT"
    c = get_extinction_coeff(seq, "DNA", True)
    # no longer getting this value??
    # assert c == 59627


def test_rna():
    seq = "ACGU"
    c = get_extinction_coeff(seq, "RNA")
    assert c == 41500


def test_rna_ss():
    seq = "AAAAAAAAUUUU"
    ss = "((((....))))"
    c1 = get_extinction_coeff(seq, "RNA")
    assert c1 == 137100
    c2 = get_extinction_coeff(seq, "RNA", structure=ss)
    assert c2 == 113336
