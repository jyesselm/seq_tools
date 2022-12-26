"""
Make sure I understand how the re module works. How to use it for sequence
and structure parsing.
"""

import re


def test_re_overlapping():
    """
    test re module udnerstanding how to do overlapping matches
    """
    # pattern = re.compile(r"GG")
    pattern = re.compile(r"(?=(GG))")
    seq = "ATCGGGATCG"
    matches = []
    for match in pattern.finditer(seq):
        matches.append([match.start(), match.end() + len(match.group(1))])
    assert matches == [[3, 5], [4, 6]]


def test_find_same_position():
    """
    test re understand how to find same position matches
    """
    seq = "AUCGGGAUCGG"
    ss = "....((...(."
    pattern_seq = re.compile(r"(?=(GG))")
    pattern_ss = re.compile(r"(?=(\(\())")
    # must save as string so sets can be compared
    matches_seq = [
        str(m.start()) + "-" + str(m.end() + len(m.group(1)))
        for m in pattern_seq.finditer(seq)
    ]
    matches_ss = [
        str(m.start()) + "-" + str(m.end() + len(m.group(1)))
        for m in pattern_ss.finditer(ss)
    ]
    inter = list(set(matches_seq).intersection(set(matches_ss)))
    # split each string in inter into a list of ints
    inter = [list(map(int, i.split("-"))) for i in inter]
    assert inter == [[4, 6]]


def test_convert_seq_and_ss_to_patterns():
    """
    test convert sequence and structure to patterns
    """
    seq_sub = r"NN"
    ss_sub = r"(("
    seq_sub = r"(?=(" + seq_sub.replace("N", r"\S") + r"))"
    ss_sub = r"(?=(" + ss_sub.replace("(", r"\(") + r"))"
    assert r"(?=(\S\S))" == seq_sub
    assert r"(?=(\(\())" == ss_sub
    #pattern_seq = re.compile(r"(?=(GG))")
    #pattern_ss = re.compile(r"(?=(\(\())")