"""
a module for calculating extinction coefficients for nucleic acids
"""

from seq_tools import dot_bracket, sequence


def get_extinction_coeff(seq, ntype, double_stranded=False, structure=None):
    """
    get the extinction coefficient for a sequence
    :param seq: sequence
    :param ntype: DNA or RNA
    :param double_stranded: is double stranded?
    :param structure: structure of the sequence in dot bracket notation
    :return: float
    """
    dna_di = {
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

    dna_mono = {"A": 15400, "C": 7400, "G": 11500, "T": 8700}

    rna_di = {
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

    rna_mono = {"A": 15400, "C": 7400, "G": 11500, "U": 9900}

    def get_mono_contribution(seq, ntype) -> float:
        """
        get the contribution of the mononucleotides to the extinction coefficient
        :param seq: sequence
        :param ntype: type of nucleic acid (DNA or RNA)
        :return: float
        """
        total = 0
        for nuc in seq[1:-1]:
            if ntype == "RNA":
                total += rna_mono[nuc]
            else:
                total += dna_mono[nuc]
        return total

    def get_di_contribution(seq, ntype) -> float:
        """
        get the contribution of the dinucleotides to the extinction coefficient
        :param seq: sequence
        :param ntype: DNA or RNA
        :return: float
        """
        total = 0
        for i in range(0, len(seq) - 1):
            distep = seq[i] + seq[i + 1]
            if ntype == "RNA":
                total += rna_di[distep]
            else:
                total += dna_di[distep]
        return total

    def get_hypochromicity_dna(seq) -> float:
        """
        get the hypochromicity of a DNA sequence
        :param seq: sequence
        :return: float
        """
        frac_at = 0
        for nuc in seq:
            if nuc in ("A", "T"):
                frac_at += 1
        frac_at /= len(seq)
        return frac_at * 0.287 + (1 - frac_at) * 0.059

    def get_hypochromicity_rna(seq, secstruct) -> float:
        """
        get the hypochromicity of an RNA sequence
        :param seq: sequence
        :param secstruct: secondary structure in dot-bracket notation
        :return: float
        """
        pairtable = dot_bracket.dotbracket_to_pairtable(secstruct)
        frac_au = 0
        frac_gc = 0
        for cur_pos in pairtable:
            if cur_pos == -1:
                continue
            pos = cur_pos
            pos2 = pairtable[pos]
            name = seq[pos] + seq[pos2]
            if name in ("AU", "UA"):
                frac_au += 1
            if name in ("GC", "CG"):
                frac_gc += 1
        frac_au /= len(seq)
        frac_gc /= len(seq)
        return frac_au * 0.26 + frac_gc * 0.059

    def get_coefficient_dna(seq, double_stranded=False) -> float:
        """
        get the extinction coefficient for a DNA sequence
        :param seq:
        :param double_stranded:
        :return: float
        """
        mono_con = get_mono_contribution(seq, "DNA")
        di_con = get_di_contribution(seq, "DNA")
        strand1 = di_con - mono_con
        if not double_stranded:
            return strand1
        rev_comp = sequence.get_reverse_complement(seq, "DNA")
        strand2 = get_di_contribution(rev_comp, "DNA") - get_mono_contribution(
            rev_comp, "DNA"
        )
        hc_val = get_hypochromicity_dna(seq)
        final = round((1 - hc_val) * (strand1 + strand2))
        return final

    def get_coefficient_rna(seq, secstruct=None):
        mono_cont = get_mono_contribution(seq, "RNA")
        di_cont = get_di_contribution(seq, "RNA")
        if secstruct is not None:
            hc_val = get_hypochromicity_rna(seq, secstruct)
            return round((1 - hc_val) * (di_cont - mono_cont))
        return di_cont - mono_cont

    if ntype == "RNA":
        return get_coefficient_rna(seq, structure)
    return get_coefficient_dna(seq, double_stranded)
