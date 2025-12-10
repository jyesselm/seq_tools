"""Core dataframe operations for sequence manipulation."""

from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import vienna


def split(df: pd.DataFrame, n_chunks: int) -> list[pd.DataFrame]:
    """Split a DataFrame into multiple chunks.

    Args:
        df: The DataFrame to be split.
        n_chunks: The number of chunks to split the DataFrame into.

    Returns:
        A list of DataFrames, each representing a chunk of the original DataFrame.
    """
    chunk_size = len(df) // n_chunks
    chunks = [df.iloc[i * chunk_size : (i + 1) * chunk_size] for i in range(n_chunks)]

    # Handle the last chunk in case the division isn't perfect
    if len(df) % n_chunks != 0:
        chunks.append(df.iloc[n_chunks * chunk_size :])

    return chunks


def run_in_parallel(
    df: pd.DataFrame, func: Callable[[pd.DataFrame], pd.DataFrame], threads: int
) -> pd.DataFrame:
    """Run a function in parallel on chunks of a DataFrame using multiple threads.

    Args:
        df: The DataFrame to process.
        func: The function to apply to each chunk of the DataFrame.
        threads: The number of threads to use for parallel processing.

    Returns:
        The combined results of applying the function to each chunk.
    """
    df_chunks = split(df, threads)
    results = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(func, chunk): chunk for chunk in df_chunks}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
    # Combine the results into a single DataFrame
    df_results = pd.concat(results)
    return df_results


def fold(df: pd.DataFrame) -> pd.DataFrame:
    """Fold each sequence in the dataframe using ViennaRNA.

    Adds columns for structure (dot-bracket notation), mfe (minimum free energy),
    and ens_defect (ensemble defect).

    Args:
        df: DataFrame containing sequences to fold.

    Returns:
        DataFrame with added 'structure', 'mfe', and 'ens_defect' columns.
    """

    def _fold(seq: str) -> pd.Series:
        v_res = vienna.fold(seq)
        return pd.Series(
            [v_res.dot_bracket, v_res.mfe, v_res.ens_defect],
            index=["structure", "mfe", "ens_defect"],
        )

    df = df.copy()
    df[["structure", "mfe", "ens_defect"]] = df["sequence"].apply(_fold)
    return df


def add(df: pd.DataFrame, p5_seq: str = "", p3_seq: str = "") -> pd.DataFrame:
    """Add a 5' and 3' sequence to the sequences in the dataframe.

    Args:
        df: DataFrame containing sequences.
        p5_seq: 5' sequence to prepend.
        p3_seq: 3' sequence to append.

    Returns:
        DataFrame with modified sequences. If structure column exists, it will be re-folded.
    """
    df = df.copy()
    df["sequence"] = df["sequence"].apply(lambda x: p5_seq + x + p3_seq)
    if "structure" in df.columns:
        df = fold(df)
    return df


def _trim_string_column(column: pd.Series, start: int, end: int) -> pd.Series:
    """Trim a string column based on start and end indices.

    Args:
        column: Series of strings to trim.
        start: Number of characters to remove from start.
        end: Number of characters to remove from end.

    Returns:
        Trimmed Series.
    """
    if start == 0 and end == 0:
        return column
    if end == 0:
        return column.str[start:]
    if start == 0:
        return column.str[:-end]
    return column.str[start:-end]


def _trim_list_column(column: pd.Series, start: int, end: int) -> pd.Series:
    """Trim a column containing lists or arrays.

    Args:
        column: Series of lists/arrays to trim.
        start: Number of elements to remove from start.
        end: Number of elements to remove from end.

    Returns:
        Trimmed Series.
    """
    if end != 0:
        return column.apply(lambda x: x[start:-end])
    return column.apply(lambda x: x[start:])


def trim(
    df: pd.DataFrame, start: int, end: int, extra_columns: list[str] | None = None
) -> pd.DataFrame:
    """Trim the sequence and other columns to the given start and end indices.

    Args:
        df: A DataFrame with 'sequence' and optionally 'structure' columns.
        start: The start index for trimming (number of bases to remove from 5' end).
        end: The end index for trimming (number of bases to remove from 3' end).
        extra_columns: Additional columns to trim along with sequence/structure.

    Returns:
        A trimmed DataFrame with columns adjusted to the specified indices.
    """
    if extra_columns is None:
        extra_columns = []

    df = df.copy()
    trim_columns = ["sequence"]
    if "structure" in df.columns:
        trim_columns.append("structure")
    trim_columns.extend(extra_columns)

    for col in trim_columns:
        if col not in df.columns:
            continue
        first_val = df.iloc[0][col]
        if isinstance(first_val, (list, pd.arrays.NumpyExtensionArray)):
            df[col] = _trim_list_column(df[col], start, end)
        else:
            df[col] = _trim_string_column(df[col], start, end)

    return df


def get_length(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate the length of each sequence in the dataframe.

    Args:
        df: DataFrame containing sequences.

    Returns:
        DataFrame with added 'length' column.
    """
    df = df.copy()
    df["length"] = df["sequence"].apply(len)
    return df


def get_default_names(df: pd.DataFrame) -> pd.DataFrame:
    """Add default names to dataframe if not already present.

    Creates names in the form 'seq_0', 'seq_1', etc. based on row indices.

    Args:
        df: DataFrame without a 'name' column.

    Returns:
        DataFrame with added 'name' column.

    Raises:
        ValueError: If dataframe already has a 'name' column.
    """
    if "name" in df.columns:
        raise ValueError("Dataframe already has names")
    df = df.copy()
    df["name"] = df.index.map(lambda x: f"seq_{x}")
    return df
