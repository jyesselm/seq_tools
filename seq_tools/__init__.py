"""
seq_tools - Tools for working with nucleic acid sequences

This package provides utilities for manipulating and analyzing nucleic acid
sequences (DNA and RNA) both as individual sequences and in batch operations
using pandas DataFrames.
"""

__author__ = "Joe Yesselman"
__email__ = "jyesselm@unl.edu"
__version__ = "0.10.0"

# Single sequence functions available via submodule
from . import extinction_coeff

# DataFrame-level functions (for batch operations)
from .dataframe import (
    add,
    fold,
    trim,
    to_fasta,
    to_opool,
    transcribe,
    determine_ntype,
    get_extinction_coeff,
    get_molecular_weight,
    get_reverse_complement,
    get_default_names,
    get_length,
    calc_edit_distance,
    calc_edit_distance_parallel,
    generate_mutated_sequences,
    generate_random_sequences,
    has_t7_promoter,
    has_5p_sequence,
    has_3p_sequence,
    has_sequence,
    to_dna,
    to_rna,
    to_dna_template,
)

# Utility functions
from .utils import (
    sequence_to_dataframe,
    sequences_to_dataframe,
    dataframe_to_sequences,
    get_resources_path,
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
    # Single sequence functions (via submodule)
    "extinction_coeff",
    # DataFrame-level functions (batch operations)
    "add",
    "fold",
    "trim",
    "to_fasta",
    "to_opool",
    "transcribe",
    "determine_ntype",
    "get_extinction_coeff",
    "get_molecular_weight",
    "get_reverse_complement",
    "get_length",
    "calc_edit_distance",
    "calc_edit_distance_parallel",
    "generate_mutated_sequences",
    "generate_random_sequences",
    "has_t7_promoter",
    "has_5p_sequence",
    "has_3p_sequence",
    "has_sequence",
    "to_dna",
    "to_rna",
    "to_dna_template",
    # Utility functions
    "sequence_to_dataframe",
    "sequences_to_dataframe",
    "dataframe_to_sequences",
    "get_resources_path",
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
