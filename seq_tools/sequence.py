"""Simple functions for gathering information about a sequence."""

from seq_tools.config import DNA_MW, RC_DNA, RC_RNA, RNA_MW, T7_PROMOTER


def get_max_stretch(seq: str) -> float:
    """Compute max stretch of the same letter in string.

    Args:
        seq: The sequence to analyze.

    Returns:
        The maximum number of consecutive identical characters.
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


def get_molecular_weight(
    seq: str, ntype: str = "DNA", double_stranded: bool = False
) -> float:
    """Calculate the molecular weight of a nucleic acid sequence.

    Args:
        seq: The nucleotide sequence.
        ntype: Type of nucleic acid: "DNA" or "RNA". Defaults to "DNA".
        double_stranded: Whether the sequence is double-stranded. Defaults to False.

    Returns:
        Molecular weight in Daltons.

    Examples:
        >>> get_molecular_weight("ATCG", "DNA")
        1331.6
        >>> get_molecular_weight("AUCG", "RNA", double_stranded=True)
        2726.8
    """

    def compute_mw(seq: str, ntype: str) -> float:
        molecular_weight = 0.0
        for nuc in seq:
            if ntype == "RNA":
                molecular_weight += RNA_MW[nuc]
            else:
                molecular_weight += DNA_MW[nuc]
        return molecular_weight

    # enforce RNA or DNA typing
    if ntype == "RNA":
        seq = to_rna(seq)
    else:
        seq = to_dna(seq)
    molecular_weight = compute_mw(seq, ntype)
    if double_stranded:
        rev_comp = get_reverse_complement(seq, ntype)
        molecular_weight += compute_mw(rev_comp, ntype)
    return molecular_weight


def get_reverse_complement(seq: str, ntype: str = "DNA") -> str:
    """Calculate the reverse complement of a nucleic acid sequence.

    Args:
        seq: Sequence to reverse complement.
        ntype: Type of nucleic acid: "DNA" or "RNA". Defaults to "DNA".

    Returns:
        Reverse complement of the sequence.

    Examples:
        >>> get_reverse_complement("ATCG", "DNA")
        'CGAT'
        >>> get_reverse_complement("AUCG", "RNA")
        'CGAU'
    """
    if ntype == "RNA":
        seq = to_rna(seq)
    else:
        seq = to_dna(seq)
    rev_comp = ""
    for nuc in seq:
        if ntype == "RNA":
            rev_comp += RC_RNA[nuc]
        else:
            rev_comp += RC_DNA[nuc]
    return rev_comp[::-1]


def to_dna(seq: str) -> str:
    """Convert RNA sequence to DNA (replaces U with T).

    Args:
        seq: RNA sequence.

    Returns:
        DNA sequence.

    Examples:
        >>> to_dna("AUCG")
        'ATCG'
    """
    return seq.replace("U", "T")


def to_dna_template(seq: str) -> str:
    """Convert RNA sequence to DNA template with T7 promoter.

    Args:
        seq: RNA sequence.

    Returns:
        DNA template sequence with T7 promoter at 5' end.

    Examples:
        >>> to_dna_template("AUCG")
        'TTCTAATACGACTCACTATAATCG'
    """
    return T7_PROMOTER + to_dna(seq)


def to_rna(seq: str) -> str:
    """Convert DNA sequence to RNA (replaces T with U).

    Args:
        seq: DNA sequence.

    Returns:
        RNA sequence.

    Examples:
        >>> to_rna("ATCG")
        'AUCG'
    """
    return seq.replace("T", "U")
