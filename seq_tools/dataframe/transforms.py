"""Sequence transformation operations for dataframes."""

import pandas as pd

from seq_tools import extinction_coeff, sequence


def to_dna(df: pd.DataFrame) -> pd.DataFrame:
    """Convert each sequence in dataframe to DNA.

    Replaces U with T in all sequences. Removes structure column if present.

    Args:
        df: DataFrame containing sequences (presumably RNA).

    Returns:
        DataFrame with sequences converted to DNA.
    """
    df = df.copy()
    df["sequence"] = df["sequence"].apply(sequence.to_dna)
    if "structure" in df.columns:
        df = df.drop(columns=["structure"])
    return df


def to_rna(df: pd.DataFrame) -> pd.DataFrame:
    """Convert each sequence in dataframe to RNA.

    Replaces T with U in all sequences.

    Args:
        df: DataFrame containing sequences (presumably DNA).

    Returns:
        DataFrame with sequences converted to RNA.
    """
    df = df.copy()
    df["sequence"] = df["sequence"].apply(sequence.to_rna)
    return df


def to_dna_template(df: pd.DataFrame) -> pd.DataFrame:
    """Convert each sequence in dataframe to DNA template.

    Replaces U with T and reverse complements sequences. Removes structure column if present.

    Args:
        df: DataFrame containing sequences.

    Returns:
        DataFrame with sequences converted to DNA template.
    """
    df = df.copy()
    df["sequence"] = df["sequence"].apply(sequence.to_dna_template)
    if "structure" in df.columns:
        df = df.drop(columns=["structure"])
    return df


def get_reverse_complement(df: pd.DataFrame, ntype: str) -> pd.DataFrame:
    """Reverse complement each sequence in the dataframe.

    Args:
        df: DataFrame containing sequences.
        ntype: Nucleotide type, "RNA" or "DNA".

    Returns:
        DataFrame with added 'rev_comp' column containing reverse complements.
    """
    df = df.copy()
    df["rev_comp"] = df["sequence"].apply(
        lambda x: sequence.get_reverse_complement(x, ntype)
    )
    return df


def get_molecular_weight(
    df: pd.DataFrame, ntype: str, double_stranded: bool = False
) -> pd.DataFrame:
    """Calculate the molecular weight for each sequence in the dataframe.

    Args:
        df: DataFrame containing sequences.
        ntype: Nucleotide type, "RNA" or "DNA".
        double_stranded: Whether sequences are double-stranded.

    Returns:
        DataFrame with added 'mw' column containing molecular weights.
    """
    df = df.copy()
    df["mw"] = df["sequence"].apply(
        lambda x: sequence.get_molecular_weight(x, ntype, double_stranded)
    )
    return df


def get_extinction_coeff(
    df: pd.DataFrame, ntype: str, double_stranded: bool = False
) -> pd.DataFrame:
    """Calculate the extinction coefficient for each sequence in the dataframe.

    For RNA sequences with structure information, uses structure-aware calculation.
    Otherwise, uses standard sequence-based calculation.

    Args:
        df: DataFrame containing sequences.
        ntype: Nucleotide type, "RNA" or "DNA".
        double_stranded: Whether sequences are double-stranded.

    Returns:
        DataFrame with added 'extinction_coeff' column.
    """

    def compute_w_struct(row: pd.Series) -> float:
        """Compute the extinction coefficient for a sequence with a structure.

        Args:
            row: DataFrame row containing 'sequence' and 'structure' columns.

        Returns:
            Extinction coefficient value.
        """
        return extinction_coeff.get_extinction_coeff(
            row["sequence"], ntype, double_stranded, row["structure"]
        )

    df = df.copy()
    if ntype == "RNA" and "structure" in df.columns:
        df["extinction_coeff"] = df.apply(compute_w_struct, axis=1)
    else:
        df["extinction_coeff"] = df["sequence"].apply(
            lambda x: extinction_coeff.get_extinction_coeff(x, ntype, double_stranded)
        )
    return df
