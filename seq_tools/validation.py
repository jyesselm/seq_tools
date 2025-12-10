"""Validation functions for sequences and dataframes."""

from typing import Optional

import pandas as pd


def validate_sequence(
    seq: str, ntype: Optional[str] = None, allow_ambiguous: bool = True
) -> None:
    """Validate sequence format.

    Args:
        seq: The sequence to validate.
        ntype: Expected nucleotide type: "DNA" or "RNA". If None, auto-detects.
        allow_ambiguous: Whether to allow ambiguous nucleotides (N). Defaults to True.

    Raises:
        ValueError: If sequence is invalid or empty.
    """
    if not seq:
        raise ValueError("Sequence cannot be empty")

    valid_dna = set("ATCGN") if allow_ambiguous else set("ATCG")
    valid_rna = set("AUCGN") if allow_ambiguous else set("AUCG")

    seq_upper = seq.upper()

    if ntype == "DNA":
        invalid = set(seq_upper) - valid_dna
        if invalid:
            raise ValueError(
                f"Invalid DNA sequence: contains invalid characters: {invalid}"
            )
    elif ntype == "RNA":
        invalid = set(seq_upper) - valid_rna
        if invalid:
            raise ValueError(
                f"Invalid RNA sequence: contains invalid characters: {invalid}"
            )
    else:
        # Auto-detect: must be either DNA or RNA
        invalid_dna = set(seq_upper) - valid_dna
        invalid_rna = set(seq_upper) - valid_rna
        if invalid_dna and invalid_rna:
            raise ValueError(
                f"Invalid sequence: contains characters that are neither DNA nor RNA: "
                f"{set(seq_upper) - valid_dna - valid_rna}"
            )


def validate_dataframe(df: pd.DataFrame, require_name: bool = True) -> None:
    """Validate a dataframe to have required columns.

    Args:
        df: DataFrame to validate.
        require_name: Whether to require a 'name' column. Defaults to True.

    Raises:
        ValueError: If required columns are missing.
    """
    if "sequence" not in df.columns:
        raise ValueError("DataFrame must have a 'sequence' column")

    if require_name and "name" not in df.columns:
        raise ValueError("DataFrame must have a 'name' column")


def ensure_name_column(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure dataframe has a 'name' column, adding default names if missing.

    Args:
        df: DataFrame to process.

    Returns:
        DataFrame with 'name' column guaranteed to exist.
    """
    if "name" not in df.columns:
        df = df.copy()
        df["name"] = [f"seq_{i}" for i in range(len(df))]
    return df
