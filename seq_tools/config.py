"""
Configuration constants for seq_tools
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
