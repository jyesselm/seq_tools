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

# Configuration constants
from .config import (
    DEFAULT_DNA_NTS,
    DEFAULT_RNA_NTS,
    DNA_MW,
    RC_DNA,
    RC_RNA,
    RNA_MW,
    T7_PROMOTER,
)

# DataFrame-level functions (for batch operations)
# Note: These are imported from the dataframe package (refactored from single module)
from .dataframe import (
    add,
    calc_edit_distance,
    calc_edit_distance_parallel,
    determine_ntype,
    fold,
    generate_mutated_sequences,
    generate_random_sequences,
    get_extinction_coeff,
    get_length,
    get_molecular_weight,
    get_reverse_complement,
    has_3p_sequence,
    has_5p_sequence,
    has_sequence,
    has_t7_promoter,
    to_dna,
    to_dna_template,
    to_fasta,
    to_opool,
    to_rna,
    transcribe,
    trim,
)

# Utility functions
from .utils import (
    dataframe_to_sequences,
    get_resources_path,
    sequence_to_dataframe,
    sequences_to_dataframe,
)

# Validation functions
from .validation import (
    ensure_name_column,
    validate_dataframe,
    validate_sequence,
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
