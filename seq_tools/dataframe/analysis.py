"""Analysis functions for sequence dataframes."""

import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Optional

import editdistance
import numpy as np
import pandas as pd


def calc_edit_distance(df: pd.DataFrame) -> float:
    """Calculate the average minimum edit distance between sequences in the dataframe.

    For each sequence, finds the minimum edit distance to any other sequence,
    then returns the average of these minimum distances.

    Args:
        df: DataFrame containing sequences.

    Returns:
        Average minimum edit distance. Returns 0.0 if dataframe has only one sequence.
    """
    if len(df) == 1:
        return 0.0
    scores = [100 for _ in range(len(df))]
    sequences = list(df["sequence"])
    for i, seq1 in enumerate(sequences):
        for j, seq2 in enumerate(sequences):
            if i >= j:
                continue
            diff = editdistance.eval(seq1, seq2)
            if scores[i] > diff:
                scores[i] = diff
            if scores[j] > diff:
                scores[j] = diff
    avg = float(np.mean(scores))
    return avg


def _compute_pairwise_edit_distances(
    chunk: list[tuple[int, int, str, str]],
) -> dict[int, int]:
    """Compute pairwise edit distances for a chunk of sequence pairs.

    Helper function used for parallel processing.

    Args:
        chunk: List of tuples (i, j, seq1, seq2) representing pairs to compare.

    Returns:
        Dictionary mapping sequence indices to their minimum edit distances.
    """
    results = {}
    for i, j, seq1, seq2 in chunk:
        diff = editdistance.eval(seq1, seq2)
        if i not in results:
            results[i] = diff
        else:
            results[i] = min(results[i], diff)
        if j not in results:
            results[j] = diff
        else:
            results[j] = min(results[j], diff)
    return results


def calc_edit_distance_parallel(
    df: pd.DataFrame, n_workers: Optional[int] = None, use_threads: bool = False
) -> float:
    """Calculate the average minimum edit distance using parallel processing.

    For each sequence, finds the minimum edit distance to any other sequence,
    then returns the average of these minimum distances. Uses parallel processing
    to speed up computation for large datasets.

    Args:
        df: DataFrame containing sequences.
        n_workers: Number of workers to use. If None, uses the number of CPU cores.
        use_threads: If True, use threads instead of processes.

    Returns:
        Average minimum edit distance. Returns 0.0 if dataframe has only one sequence.
    """
    if len(df) == 1:
        return 0.0

    sequences = list(df["sequence"])
    n = len(sequences)

    # Generate all pairs to compare
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((i, j, sequences[i], sequences[j]))

    # Split pairs into chunks for parallel processing
    if n_workers is None:
        n_workers = os.cpu_count() or 1

    chunk_size = max(1, len(pairs) // n_workers)
    chunks = [pairs[i : i + chunk_size] for i in range(0, len(pairs), chunk_size)]

    # Process chunks in parallel
    scores = [100 for _ in range(n)]
    executor_class = ThreadPoolExecutor if use_threads else ProcessPoolExecutor
    with executor_class(max_workers=n_workers) as executor:
        futures = [
            executor.submit(_compute_pairwise_edit_distances, chunk) for chunk in chunks
        ]
        for future in as_completed(futures):
            chunk_results = future.result()
            # Update scores with minimum values
            for idx, dist in chunk_results.items():
                scores[idx] = min(scores[idx], dist)

    avg = float(np.mean(scores))
    return avg


def determine_ntype(df: pd.DataFrame) -> str:
    """Determine the nucleotide type of the sequences in the dataframe.

    Analyzes sequences for presence of T (DNA) or U (RNA) nucleotides.
    Raises an error if both types are found in longer sequences.

    Args:
        df: DataFrame containing sequences.

    Returns:
        "RNA" or "DNA" based on nucleotide content.

    Raises:
        ValueError: If both DNA and RNA sequences are detected in longer sequences.
    """
    results = []
    for _, row in df.iterrows():
        ntype = "UNCERTAIN"
        if row["sequence"].count("T") > 0:
            ntype = "DNA"
        elif row["sequence"].count("U") > 0:
            ntype = "RNA"
        results.append(ntype)
    if df["sequence"].str.len().mean() > 10:
        if results.count("DNA") > 0 and results.count("RNA") > 0:
            raise ValueError("Cannot determine nucleotide type")
    if results.count("RNA") > 0:
        return "RNA"
    return "DNA"
