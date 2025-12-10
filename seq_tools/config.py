"""Configuration constants for seq_tools.

This module defines nucleotide types, molecular weights, reverse complement
mappings, and other constants used throughout seq_tools.

Attributes:
    T7_PROMOTER: Standard T7 RNA polymerase promoter sequence.
    DEFAULT_DNA_NTS: List of standard DNA nucleotides.
    DEFAULT_RNA_NTS: List of standard RNA nucleotides.
    VALID_DNA_NTS: Set of valid DNA nucleotides including ambiguous bases.
    VALID_RNA_NTS: Set of valid RNA nucleotides including ambiguous bases.
    RNA_MW: Dictionary mapping RNA nucleotides to molecular weights in Daltons.
    DNA_MW: Dictionary mapping DNA nucleotides to molecular weights in Daltons.
    RC_DNA: Dictionary for DNA reverse complement mapping.
    RC_RNA: Dictionary for RNA reverse complement mapping.
"""

# T7 promoter sequence
T7_PROMOTER = "TTCTAATACGACTCACTATA"

# Default nucleotide types
DEFAULT_DNA_NTS = ["A", "T", "G", "C"]
DEFAULT_RNA_NTS = ["A", "U", "G", "C"]

# Valid nucleotide sets
VALID_DNA_NTS = set("ATCGN")
VALID_RNA_NTS = set("AUCGN")

# Molecular weights (RNA) - in Daltons
RNA_MW = {"A": 347.2, "C": 323.2, "G": 363.2, "U": 324.2}

# Molecular weights (DNA) - in Daltons
DNA_MW = {"A": 331.2, "C": 307.2, "G": 347.2, "T": 322.2}

# Reverse complement mappings
RC_DNA = {"A": "T", "T": "A", "G": "C", "C": "G"}
RC_RNA = {"A": "U", "U": "A", "G": "C", "C": "G"}
