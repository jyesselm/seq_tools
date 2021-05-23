from seq_tools import dot_bracket, sequence
import collections as col

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


def get_mono_contribution(seq, type):
    total = 0
    for e in seq[1:-1]:
        if type == "RNA":
            total += rna_mono[e]
        else:
            total += dna_mono[e]
    return total


def get_di_contribution(seq, type):
    total = 0
    for i in range(0, len(seq) - 1):
        distep = seq[i] + seq[i + 1]
        if type == "RNA":
            total += rna_di[distep]
        else:
            total += dna_di[distep]
    return total


def get_hypochromicity_dna(seq):
    frac_at = 0
    for e in seq:
        if e == "A" or e == "T":
            frac_at += 1
    frac_at /= len(seq)
    return frac_at * 0.287 + (1 - frac_at) * 0.059


def get_hypochromicity_rna(seq, ss):
    pairtable = dot_bracket.dotbracket_to_pairtable(ss)
    frac_au = 0
    frac_gc = 0
    for i in range(0, len(pairtable)):
        if pairtable[i] == -1:
            continue
        pos = pairtable[i]
        pos2 = pairtable[pos]
        name = seq[pos] + seq[pos2]
        if name == "AU" or name == "UA":
            frac_au += 1
        if name == "GC" or name == "CG":
            frac_gc += 1
    frac_au /= len(seq)
    frac_gc /= len(seq)
    return frac_au * 0.26 + frac_gc * 0.059


def get_coefficient_dna(seq, ds=False):
    mono = get_mono_contribution(seq, "DNA")
    di = get_di_contribution(seq, "DNA")
    strand1 = di - mono
    if not ds:
        return strand1
    rc = sequence.get_reverse_complement(seq, "DNA")
    strand2 = get_di_contribution(rc, "DNA") - get_mono_contribution(rc, "DNA")
    hc = get_hypochromicity_dna(seq)
    final = round((1 - hc) * (strand1 + strand2))
    return final


def get_coefficient_rna(seq, ss=None):
    mono = get_mono_contribution(seq, "RNA")
    di = get_di_contribution(seq, "RNA")
    if ss is not None:
        hc = get_hypochromicity_rna(seq, ss)
        return round((1 - hc) * (di - mono))
    else:
        return di - mono
