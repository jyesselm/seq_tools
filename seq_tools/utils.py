"""Utility functions for converting between sequences and DataFrames."""

from pathlib import Path
from typing import Optional

import pandas as pd

try:
    # Python 3.9+
    from importlib.resources import files

    _HAS_FILES = True
except ImportError:
    try:
        # Python 3.7-3.8 with importlib_resources backport
        from importlib_resources import files

        _HAS_FILES = True
    except ImportError:
        # Fallback for older Python or if importlib_resources not available

        _HAS_FILES = False


def get_resources_path() -> Path:
    """Get the path to the seq_tools resources directory.

    This function provides a reliable way to access package resource files
    whether the package is installed or run in development mode.

    Returns:
        Path object pointing to the seq_tools/resources directory.

    Examples:
        >>> resource_path = get_resources_path()
        >>> p5_file = resource_path / "p5_sequences.csv"
        >>> df = pd.read_csv(p5_file)
    """
    if _HAS_FILES:
        return Path(files("seq_tools")) / "resources"
    else:
        # Fallback: find resources relative to this file
        import seq_tools

        seq_tools_path = Path(seq_tools.__file__).parent
        return seq_tools_path / "resources"


def sequence_to_dataframe(seq: str, name: str = "sequence") -> pd.DataFrame:
    """Convert a single sequence to a DataFrame.

    Args:
        seq: The nucleotide sequence.
        name: Name for the sequence. Defaults to "sequence".

    Returns:
        DataFrame with 'name' and 'sequence' columns.

    Examples:
        >>> df = sequence_to_dataframe("ATCG", name="seq1")
        >>> len(df)
        1
        >>> df.iloc[0]['sequence']
        'ATCG'
    """
    return pd.DataFrame({"name": [name], "sequence": [seq]})


def sequences_to_dataframe(
    sequences: list[str], names: Optional[list[str]] = None
) -> pd.DataFrame:
    """Convert multiple sequences to a DataFrame.

    Args:
        sequences: List of nucleotide sequences.
        names: List of names for sequences. If None, generates default names.

    Returns:
        DataFrame with 'name' and 'sequence' columns.

    Raises:
        ValueError: If number of sequences and names don't match.

    Examples:
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


def dataframe_to_sequences(df: pd.DataFrame) -> list[dict[str, str]]:
    """Convert DataFrame to list of sequence dictionaries.

    Args:
        df: DataFrame with 'name' and 'sequence' columns.

    Returns:
        List of dictionaries with 'name' and 'sequence' keys.

    Examples:
        >>> df = pd.DataFrame({"name": ["seq1"], "sequence": ["ATCG"]})
        >>> result = dataframe_to_sequences(df)
        >>> result[0]['sequence']
        'ATCG'
    """
    return df.to_dict("records")
