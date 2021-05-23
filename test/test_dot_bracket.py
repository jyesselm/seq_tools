import pytest
from seq_tools import dot_bracket


def test_dot_bracket():
    ss = "((()))"
    pair_table = dot_bracket.dotbracket_to_pairtable(ss)
    assert pair_table[0] == 5
    assert pair_table[-1] == 0


def test_dot_bracket_2():
    ss = "((((....))))"
    pair_table = dot_bracket.dotbracket_to_pairtable(ss)
    assert pair_table[4] == -1
