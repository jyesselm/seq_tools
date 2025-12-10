"""Package for working with dataframes that contain nucleotide sequences."""

# Analysis functions
from seq_tools.dataframe.analysis import (
    calc_edit_distance,
    calc_edit_distance_parallel,
    determine_ntype,
)

# Core operations
from seq_tools.dataframe.core import (
    add,
    fold,
    get_default_names,
    get_length,
    run_in_parallel,
    split,
    trim,
)

# Sequence generators
from seq_tools.dataframe.generators import (
    generate_mutated_sequences,
    generate_random_sequences,
)

# I/O operations
from seq_tools.dataframe.io import to_fasta, to_opool

# Primer and common sequence operations
from seq_tools.dataframe.primers import (
    add_common_seqs,
    find_longest_common_prefix,
    find_longest_common_suffix,
    has_3p_sequence,
    has_5p_sequence,
    has_sequence,
    has_t7_promoter,
    remove_common_p5_p3,
    remove_common_p5_p3_by_structure,
    remove_common_seqs,
    transcribe,
    trim_p5_and_p3,
)

# Transformation operations
from seq_tools.dataframe.transforms import (
    get_extinction_coeff,
    get_molecular_weight,
    get_reverse_complement,
    to_dna,
    to_dna_template,
    to_rna,
)

__all__ = [
    # Core operations
    "add",
    "fold",
    "trim",
    "get_length",
    "get_default_names",
    "split",
    "run_in_parallel",
    # Transforms
    "to_dna",
    "to_rna",
    "to_dna_template",
    "get_reverse_complement",
    "get_molecular_weight",
    "get_extinction_coeff",
    # Analysis
    "calc_edit_distance",
    "calc_edit_distance_parallel",
    "determine_ntype",
    # I/O
    "to_fasta",
    "to_opool",
    # Primers
    "has_5p_sequence",
    "has_3p_sequence",
    "has_sequence",
    "has_t7_promoter",
    "find_longest_common_prefix",
    "find_longest_common_suffix",
    "trim_p5_and_p3",
    "remove_common_p5_p3",
    "remove_common_p5_p3_by_structure",
    "remove_common_seqs",
    "add_common_seqs",
    "transcribe",
    # Generators
    "generate_mutated_sequences",
    "generate_random_sequences",
]
