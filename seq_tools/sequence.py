"""
simple functions for gathering information about a sequence.
"""


def to_dna(seq) -> str:
    """
    Convert RNA sequence to DNA
    :param seq: RNA sequence
    :return: DNA sequence
    """
    return seq.replace("U", "T")


def to_rna(seq) -> str:
    """
    Convert DNA sequence to RNA
    :param seq: DNA sequence
    :return: RNA sequence
    """
    return seq.replace("T", "U")


def get_reverse_complement(seq, ntype="DNA") -> str:
    """
    returns the reverse complement of a sequence
    :param seq: sequence to reverse complement
    :param ntype: type of sequence (DNA or RNA)
    :return: reverse complement of sequence
    """
    if ntype == "RNA":
        seq = to_rna(seq)
    else:
        seq = to_dna(seq)
    rev_comp = ""
    rc_dna = {"A": "T", "T": "A", "G": "C", "C": "G"}
    rc_rna = {"A": "U", "U": "A", "G": "C", "C": "G"}
    for nuc in seq:
        if ntype == "RNA":
            rev_comp += rc_rna[nuc]
        else:
            rev_comp += rc_dna[nuc]
    return rev_comp[::-1]


mw_rna = {"A": 347.2, "C": 323.2, "G": 363.2, "U": 324.2}
mw_dna = {"A": 331.2, "C": 307.2, "G": 347.2, "T": 322.2}


def get_molecular_weight(seq, type="DNA", double_stranded=False):
    total = 0
    for e in seq:
        if type == "RNA":
            total += mw_rna[e]
        else:
            total += mw_dna[e]
    if double_stranded:
        rc = get_reverse_complement(seq, type)
        for e in rc:
            if type == "RNA":
                total += mw_rna[e]
            else:
                total += mw_dna[e]
    return total


def get_max_stretch(seq):
    """returns max stretch of the same letter in string"""
    max_stretch = 0
    last = None
    for n in seq:
        if last == None:
            last = n
            stretch = 0
            continue
        if n == last:
            stretch += 1
            if stretch > max_stretch:
                max_stretch = stretch
        else:
            stretch = 0
        last = n
    return max_stretch
