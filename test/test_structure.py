"""
test structure module for seq_tools
"""
import pytest
from seq_tools.structure import SequenceStructure, find


def test_init():
    """
    test initation of SequenceStructure object
    """
    ss = SequenceStructure("ATCG", "....")
    assert ss.sequence == "ATCG"
    assert ss.structure == "...."
    # test that the sequence and structure are the same length
    with pytest.raises(ValueError):
        SequenceStructure("ATCG", "....(")


def test_split_strands():
    """
    test that split_strands returns a list of SequenceStructure objects
    """
    ss = SequenceStructure("ATCG&ATCG", "....&....")
    assert len(ss.split_strands()) == 2
    assert isinstance(ss.split_strands()[0], SequenceStructure)


def test_join():
    """
    test that join returns a SequenceStructure object
    """
    ss1 = SequenceStructure("ATCG", "....")
    ss2 = SequenceStructure("ATCG", "....")
    ss = ss1.join(ss2)
    assert isinstance(ss, SequenceStructure)
    assert ss.sequence == "ATCG&ATCG"
    assert ss.structure == "....&...."


def test_find():
    """
    test that find returns the correct index
    """
    struct = SequenceStructure("GGGAAACCC", "(((...)))")
    sub = SequenceStructure("GAAAC", "(...)")
    r = find(struct, sub)
    assert r == [([2, 7],)]


def test_find_2():
    """
    test that find returns the correct index
    """
    struct = SequenceStructure("GGGAAACCC", "(((...)))")
    sub = SequenceStructure("GGG&CCC", "(((&)))")
    r = find(struct, sub)
    assert r == [([0, 3], [6, 9])]
    assert len(r) == 1
    # should have multiple solutions
    sub = SequenceStructure("GG&CC", "((&))")
    r = find(struct, sub)
    assert len(r) == 4


def test_real_solution():
    """
    test that find returns the correct index
    """
    seq = (
        "GGGCUUCGGCCCACUGUCUAACGAGGAAACUUUGUUCAGAUGGAUAUUUCGUCAAUCUCGAGUAGGGA"
        "UUGAUAGAAAUAAGCGUGACGCGUCAGAGAAACUCUGACCAUCAUGUACAAAGAAACAACAACAACAAC"
    )
    ss = (
        "((((....)))).(((((((((((((...))))))).)))))).(((((((((((((((.....))))))"
        "))).)))))).((((((...(((((((...)))))))..))))))......................"
    )
    struct = SequenceStructure(seq, ss)
    sub = SequenceStructure("UCUAAC&GUUCAGA", "((((((&))).)))")
    r = find(struct, sub)
    assert r == [([16, 22], [33, 40])]
    bns = r[0]
    bmin, bmax = bns[0]
    assert struct[bmin:bmax].sequence == "UCUAAC"
    assert struct[bmin:bmax].structure == "(((((("
