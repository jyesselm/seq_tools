"""A module for calculating extinction coefficients for nucleic acids."""

from seq_tools import dot_bracket, sequence

# DNA dinucleotide extinction coefficients
DNA_DINUCLEOTIDE = {
    "AA": 27400,
    "AC": 21200,
    "AG": 25000,
    "AT": 22800,
    "CA": 21200,
    "CC": 14600,
    "CG": 18000,
    "CT": 15200,
    "GA": 25200,
    "GC": 17600,
    "GG": 21600,
    "GT": 20000,
    "TA": 23400,
    "TC": 16200,
    "TG": 19000,
    "TT": 16800,
}

# DNA mononucleotide extinction coefficients
DNA_MONONUCLEOTIDE = {"A": 15400, "C": 7400, "G": 11500, "T": 8700}

# RNA dinucleotide extinction coefficients
RNA_DINUCLEOTIDE = {
    "AA": 27400,
    "AC": 21200,
    "AG": 25000,
    "AU": 24000,
    "CA": 21200,
    "CC": 14600,
    "CG": 18000,
    "CU": 16200,
    "GA": 25200,
    "GC": 17600,
    "GG": 21600,
    "GU": 21200,
    "UA": 24600,
    "UC": 17200,
    "UG": 20000,
    "UU": 19600,
}

# RNA mononucleotide extinction coefficients
RNA_MONONUCLEOTIDE = {"A": 15400, "C": 7400, "G": 11500, "U": 9900}


def _get_mono_contribution(seq: str, ntype: str) -> float:
    """Calculate the contribution of mononucleotides to extinction coefficient.

    Args:
        seq: Nucleotide sequence.
        ntype: Nucleic acid type ("DNA" or "RNA").

    Returns:
        Sum of mononucleotide contributions.
    """
    mono_dict = RNA_MONONUCLEOTIDE if ntype == "RNA" else DNA_MONONUCLEOTIDE
    return sum(mono_dict[nuc] for nuc in seq[1:-1])


def _get_di_contribution(seq: str, ntype: str) -> float:
    """Calculate the contribution of dinucleotides to extinction coefficient.

    Args:
        seq: Nucleotide sequence.
        ntype: Nucleic acid type ("DNA" or "RNA").

    Returns:
        Sum of dinucleotide contributions.
    """
    di_dict = RNA_DINUCLEOTIDE if ntype == "RNA" else DNA_DINUCLEOTIDE
    total = 0
    for i in range(len(seq) - 1):
        distep = seq[i] + seq[i + 1]
        total += di_dict[distep]
    return total


def _get_hypochromicity_dna(seq: str) -> float:
    """Calculate hypochromicity of a DNA sequence.

    Args:
        seq: DNA sequence.

    Returns:
        Hypochromicity value.
    """
    frac_at = sum(1 for nuc in seq if nuc in ("A", "T")) / len(seq)
    return frac_at * 0.287 + (1 - frac_at) * 0.059


def _count_base_pairs(seq: str, pairtable: list) -> tuple[int, int]:
    """Count AU and GC base pairs in an RNA sequence.

    Args:
        seq: RNA sequence.
        pairtable: Pair table from dot-bracket notation.

    Returns:
        Tuple of (AU count, GC count).
    """
    au_count = 0
    gc_count = 0
    for pos in pairtable:
        if pos == -1:
            continue
        pair_name = seq[pos] + seq[pairtable[pos]]
        au_count += pair_name in ("AU", "UA")
        gc_count += pair_name in ("GC", "CG")
    return au_count, gc_count


def _get_hypochromicity_rna(seq: str, secstruct: str) -> float:
    """Calculate hypochromicity of an RNA sequence with secondary structure.

    Args:
        seq: RNA sequence.
        secstruct: Secondary structure in dot-bracket notation.

    Returns:
        Hypochromicity value.
    """
    pairtable = dot_bracket.dotbracket_to_pairtable(secstruct)
    au_count, gc_count = _count_base_pairs(seq, pairtable)
    frac_au = au_count / len(seq)
    frac_gc = gc_count / len(seq)
    return frac_au * 0.26 + frac_gc * 0.059


def _calculate_dna_single_strand(seq: str) -> float:
    """Calculate extinction coefficient for single-stranded DNA.

    Args:
        seq: DNA sequence.

    Returns:
        Extinction coefficient.
    """
    mono_con = _get_mono_contribution(seq, "DNA")
    di_con = _get_di_contribution(seq, "DNA")
    return di_con - mono_con


def _calculate_dna_double_strand(seq: str) -> float:
    """Calculate extinction coefficient for double-stranded DNA.

    Args:
        seq: DNA sequence.

    Returns:
        Extinction coefficient.
    """
    strand1 = _calculate_dna_single_strand(seq)
    rev_comp = sequence.get_reverse_complement(seq, "DNA")
    strand2 = _calculate_dna_single_strand(rev_comp)
    hc_val = _get_hypochromicity_dna(seq)
    return round((1 - hc_val) * (strand1 + strand2))


def _calculate_dna_coefficient(seq: str, double_stranded: bool = False) -> float:
    """Calculate extinction coefficient for DNA sequence.

    Args:
        seq: DNA sequence.
        double_stranded: Whether the DNA is double-stranded.

    Returns:
        Extinction coefficient.
    """
    if double_stranded:
        return _calculate_dna_double_strand(seq)
    return _calculate_dna_single_strand(seq)


def _calculate_rna_coefficient(seq: str, secstruct: str = None) -> float:
    """Calculate extinction coefficient for RNA sequence.

    Args:
        seq: RNA sequence.
        secstruct: Secondary structure in dot-bracket notation. Defaults to None.

    Returns:
        Extinction coefficient.
    """
    mono_cont = _get_mono_contribution(seq, "RNA")
    di_cont = _get_di_contribution(seq, "RNA")
    base_coeff = di_cont - mono_cont

    if secstruct is not None:
        hc_val = _get_hypochromicity_rna(seq, secstruct)
        return round((1 - hc_val) * base_coeff)
    return base_coeff


def get_extinction_coeff(
    seq: str, ntype: str, double_stranded: bool = False, structure: str = None
) -> float:
    """Calculate extinction coefficient for a nucleic acid sequence.

    Args:
        seq: Nucleotide sequence.
        ntype: Nucleic acid type ("DNA" or "RNA").
        double_stranded: Whether DNA is double-stranded (ignored for RNA). Defaults to False.
        structure: Secondary structure in dot-bracket notation (for RNA only). Defaults to None.

    Returns:
        Extinction coefficient.
    """
    if ntype == "RNA":
        return _calculate_rna_coefficient(seq, structure)
    return _calculate_dna_coefficient(seq, double_stranded)
