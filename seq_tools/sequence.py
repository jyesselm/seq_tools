"""
simple functions for gathering information about a sequence.
"""


def get_max_stretch(seq) -> float:
    """
    computes max stretch of the same letter in string
    """
    max_stretch = 0
    current_stretch = 0
    for i, nuc in enumerate(seq):
        if i == 0:
            current_stretch += 1
        else:
            if nuc == seq[i - 1]:
                current_stretch += 1
            else:
                if current_stretch > max_stretch:
                    max_stretch = current_stretch
                current_stretch = 1
    if current_stretch > max_stretch:
        max_stretch = current_stretch
    return max_stretch


def get_molecular_weight(seq, ntype="DNA", double_stranded=False) -> float:
    """
    returns the molecular weight of a sequence
    :param seq: the sequence
    :param ntype: type of sequence (DNA or RNA)
    :param double_stranded: is the sequence double stranded?
    :return: float
    """
    rna_mw = {"A": 347.2, "C": 323.2, "G": 363.2, "U": 324.2}
    dna_mw = {"A": 331.2, "C": 307.2, "G": 347.2, "T": 322.2}

    def compute_mw(seq, ntype):
        molecular_weight = 0
        for nuc in seq:
            if ntype == "RNA":
                molecular_weight += rna_mw[nuc]
            else:
                molecular_weight += dna_mw[nuc]
        return molecular_weight

    # enforce RNA or DNA typing
    if ntype == "RNA":
        seq = to_rna(seq)
    else:
        seq = to_dna(seq)
    molecular_weight = compute_mw(seq, ntype)
    if double_stranded:
        rev_comp = get_reverse_complement(seq, type)
        molecular_weight += compute_mw(rev_comp, ntype)
    return molecular_weight


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


def to_dna(seq) -> str:
    """
    Convert RNA sequence to DNA
    :param seq: RNA sequence
    :return: DNA sequence
    """
    return seq.replace("U", "T")


def to_dna_template(seq) -> str:
    """
    Convert RNA sequence to DNA
    :param seq: RNA sequence
    :return: DNA sequence
    """
    return "TTCTAATACGACTCACTATA" + to_dna(seq)


def to_rna(seq) -> str:
    """
    Convert DNA sequence to RNA
    :param seq: DNA sequence
    :return: RNA sequence
    """
    return seq.replace("T", "U")
