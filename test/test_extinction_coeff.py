import pytest
from seq_tools import extinction_coeff


def test_di_dna():
    seq = "AT"
    c = extinction_coeff.get_coefficient_dna(seq)
    assert c == 22800

def test_dna():
    seq = "AAT"
    c = extinction_coeff.get_coefficient_dna(seq)
    assert c == 34800

def test_ds_dna():
    seq = "ACGT"
    c = extinction_coeff.get_coefficient_dna(seq, True)
    assert c == 59627

def test_rna():
    seq = "ACGU"
    c = extinction_coeff.get_coefficient_rna(seq)
    assert c == 41500

def test_rna_ss():
    seq = "AAAAAAAAUUUU"
    ss  = "((((....))))"
    c1 = extinction_coeff.get_coefficient_rna(seq)
    assert c1 == 137100
    c2 = extinction_coeff.get_coefficient_rna(seq, ss)
    assert c2 == 113336
