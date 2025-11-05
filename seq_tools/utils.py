"""
Utility functions for converting between sequences and DataFrames
"""

from typing import List, Dict, Optional
import pandas as pd


def sequence_to_dataframe(seq: str, name: str = "sequence") -> pd.DataFrame:
    """
    Convert a single sequence to a DataFrame.

    Parameters
    ----------
    seq : str
        The nucleotide sequence.
    name : str, optional
        Name for the sequence (default: "sequence").

    Returns
    -------
    pd.DataFrame
        DataFrame with 'name' and 'sequence' columns.

    Examples
    --------
    >>> df = sequence_to_dataframe("ATCG", name="seq1")
    >>> len(df)
    1
    >>> df.iloc[0]['sequence']
    'ATCG'
    """
    return pd.DataFrame({"name": [name], "sequence": [seq]})


def sequences_to_dataframe(
    sequences: List[str], names: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Convert multiple sequences to a DataFrame.

    Parameters
    ----------
    sequences : List[str]
        List of nucleotide sequences.
    names : List[str], optional
        List of names for sequences. If None, generates default names.

    Returns
    -------
    pd.DataFrame
        DataFrame with 'name' and 'sequence' columns.

    Examples
    --------
    >>> df = sequences_to_dataframe(["ATCG", "GCTA"])
    >>> len(df)
    2
    >>> list(df['name'])
    ['seq_0', 'seq_1']
    """
    if names is None:
        names = [f"seq_{i}" for i in range(len(sequences))]
    if len(sequences) != len(names):
        raise ValueError(
            f"Number of sequences ({len(sequences)}) must match "
            f"number of names ({len(names)})"
        )
    return pd.DataFrame({"name": names, "sequence": sequences})


def dataframe_to_sequences(df: pd.DataFrame) -> List[Dict[str, str]]:
    """
    Convert DataFrame to list of sequence dictionaries.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'name' and 'sequence' columns.

    Returns
    -------
    List[Dict[str, str]]
        List of dictionaries with 'name' and 'sequence' keys.

    Examples
    --------
    >>> df = pd.DataFrame({"name": ["seq1"], "sequence": ["ATCG"]})
    >>> result = dataframe_to_sequences(df)
    >>> result[0]['sequence']
    'ATCG'
    """
    return df.to_dict("records")

