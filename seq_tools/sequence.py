"""
simple functions for gathering information about a sequence.
"""


def get_nucleic_acid_type(seq) -> str:
    """
    Returns the type of nucleic acid in the sequence
    :param seq: sequence of dna or rna
    :return: str "DNA" or "RNA"
    """
    for nuc in seq:
        if nuc == "U":
            return "RNA"
    return "DNA"


def convert_to_rna(seq):
    new_seq = ""
    for e in seq:
        if e == "T":
            new_seq += "U"
        else:
            new_seq += e
    return new_seq


def to_dna(seq):
    """
    Convert RNA sequence to DNA
    """
    return seq.replace("U", "T")


def get_reverse_complement(seq, type="DNA"):
    complement = ""
    for e in seq:
        if e == "A" and type == "RNA":
            complement += "U"
        elif e == "A" and type == "DNA":
            complement += "T"
        elif e == "C":
            complement += "G"
        elif e == "G":
            complement += "C"
        elif e == "U":
            complement += "A"
        elif e == "T":
            complement += "A"
    return complement[::-1]


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
