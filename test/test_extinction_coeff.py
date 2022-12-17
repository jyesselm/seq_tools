"""
test extinction coefficient module
"""

from seq_tools.extinction_coeff import get_extinction_coeff


def test_ds_dna():
    """
    test double stranded DNA extinction coefficient
    """
    seq = "ACGT"
    c = get_extinction_coeff(seq, "DNA", True)
    assert c == 66656


def test_rna():
    """
    test RNA extinction coefficient
    """
    seq = "ACGU"
    c = get_extinction_coeff(seq, "RNA")
    assert c == 41500


def test_rna_ss():
    """
    test single stranded RNA extinction coefficient
    """
    seq = "AAAAAAAAUUUU"
    ss = "((((....))))"
    c = get_extinction_coeff(seq, "RNA")
    assert c == 137100
    c = get_extinction_coeff(seq, "RNA", structure=ss)
    assert c == 113336
