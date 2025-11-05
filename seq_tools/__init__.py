"""
seq_tools - Tools for working with nucleic acid sequences

This package provides utilities for manipulating and analyzing nucleic acid
sequences (DNA and RNA) both as individual sequences and in batch operations
using pandas DataFrames.
"""

__author__ = "Joe Yesselman"
__email__ = "jyesselm@unl.edu"
__version__ = "0.7.2"

# Core sequence-level functions (for single sequences)
from .sequence import (
    to_dna,
    to_rna,
    to_dna_template,
    get_reverse_complement,
    get_molecular_weight,
    get_max_stretch,
)

# Extinction coefficient functions
from .extinction_coeff import get_extinction_coeff

# Structure classes and functions
from .structure import SequenceStructure, Match, RNASegment, find, find_seq_struct

# DataFrame-level functions (for batch operations)
from .dataframe import (
    add,
    fold,
    trim,
    to_fasta,
    to_opool,
    transcribe,
    determine_ntype,
    get_extinction_coeff as get_extinction_coeff_df,
    get_molecular_weight as get_molecular_weight_df,
    get_reverse_complement as get_reverse_complement_df,
    get_length,
    calc_edit_distance,
    calc_edit_distance_parallel,
    generate_mutated_sequences,
    generate_random_sequences,
    has_t7_promoter,
    has_5p_sequence,
    has_3p_sequence,
    has_sequence,
    has_seq_struct,
    to_dna as to_dna_df,
    to_rna as to_rna_df,
    to_dna_template as to_dna_template_df,
)

# Utility functions
from .utils import (
    sequence_to_dataframe,
    sequences_to_dataframe,
    dataframe_to_sequences,
)

# Validation functions
from .validation import (
    validate_sequence,
    validate_dataframe,
    ensure_name_column,
)

# Configuration constants
from .config import (
    T7_PROMOTER,
    DEFAULT_DNA_NTS,
    DEFAULT_RNA_NTS,
    DNA_MW,
    RNA_MW,
    RC_DNA,
    RC_RNA,
)

# Explicit public API
__all__ = [
    # Version info
    "__version__",
    "__author__",
    "__email__",
    # Sequence-level functions (single sequences)
    "to_dna",
    "to_rna",
    "to_dna_template",
    "get_reverse_complement",
    "get_molecular_weight",
    "get_max_stretch",
    "get_extinction_coeff",
    # Structure
    "SequenceStructure",
    "Match",
    "RNASegment",
    "find",
    "find_seq_struct",
    # DataFrame-level functions (batch operations)
    "add",
    "fold",
    "trim",
    "to_fasta",
    "to_opool",
    "transcribe",
    "determine_ntype",
    "get_extinction_coeff_df",
    "get_molecular_weight_df",
    "get_reverse_complement_df",
    "get_length",
    "calc_edit_distance",
    "calc_edit_distance_parallel",
    "generate_mutated_sequences",
    "generate_random_sequences",
    "has_t7_promoter",
    "has_5p_sequence",
    "has_3p_sequence",
    "has_sequence",
    "has_seq_struct",
    "to_dna_df",
    "to_rna_df",
    "to_dna_template_df",
    # Utility functions
    "sequence_to_dataframe",
    "sequences_to_dataframe",
    "dataframe_to_sequences",
    # Validation functions
    "validate_sequence",
    "validate_dataframe",
    "ensure_name_column",
    # Configuration
    "T7_PROMOTER",
    "DEFAULT_DNA_NTS",
    "DEFAULT_RNA_NTS",
    "DNA_MW",
    "RNA_MW",
    "RC_DNA",
    "RC_RNA",
]
