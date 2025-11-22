"""
module for working with dataframes that contain nucleotide sequences
"""

import os
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import random
from typing import Optional

import editdistance
import numpy as np
import pandas as pd
import vienna

from seq_tools import sequence, extinction_coeff
from seq_tools.config import T7_PROMOTER, DEFAULT_DNA_NTS, DEFAULT_RNA_NTS
from seq_tools.validation import validate_dataframe, ensure_name_column
from seq_tools.utils import get_resources_path


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
        return 0
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
        return 0

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


def fold(df: pd.DataFrame) -> pd.DataFrame:
    """Fold each sequence in the dataframe using ViennaRNA.

    Adds columns for structure (dot-bracket notation), mfe (minimum free energy),
    and ens_defect (ensemble defect).

    Args:
        df: DataFrame containing sequences to fold.

    Returns:
        DataFrame with added 'structure', 'mfe', and 'ens_defect' columns.
    """

    def _fold(seq):
        v_res = vienna.fold(seq)
        return pd.Series(
            [v_res.dot_bracket, v_res.mfe, v_res.ens_defect],
            index=["structure", "mfe", "ens_defect"],
        )

    df = df.copy()
    df[["structure", "mfe", "ens_defect"]] = df["sequence"].apply(_fold)
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
    # add `name` column to dataframe in the form of `seq_1`, `seq_2`, etc.
    df = df.copy()
    df["name"] = df.index.map(lambda x: "seq_" + str(x))
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


def has_5p_sequence(df: pd.DataFrame, p5_seq: str) -> bool:
    """Check if all sequences start with the given 5' sequence.

    Args:
        df: DataFrame containing sequences.
        p5_seq: 5' sequence to check for.

    Returns:
        True if all sequences start with p5_seq, False otherwise.
    """
    return df["sequence"].str.startswith(p5_seq).all()


def has_3p_sequence(df: pd.DataFrame, p3_seq: str) -> bool:
    """Check if all sequences end with the given 3' sequence.

    Args:
        df: DataFrame containing sequences.
        p3_seq: 3' sequence to check for.

    Returns:
        True if all sequences end with p3_seq, False otherwise.
    """
    return df["sequence"].str.endswith(p3_seq).all()


def has_sequence(df: pd.DataFrame, seq: str) -> bool:
    """Check if all sequences contain the given subsequence.

    Args:
        df: DataFrame containing sequences.
        seq: Subsequence to search for.

    Returns:
        True if all sequences contain seq, False otherwise.
    """
    return df["sequence"].str.contains(seq).all()


def has_t7_promoter(df: pd.DataFrame) -> bool:
    """Check if all sequences in the dataframe have a T7 promoter.

    Args:
        df: DataFrame with sequences.

    Returns:
        True if all sequences start with T7 promoter, False otherwise.
    """
    has_t7 = df[df["sequence"].str.startswith(T7_PROMOTER)]
    if len(has_t7) != len(df):
        return False
    return True


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


def to_fasta(df: pd.DataFrame, filename: str) -> None:
    """Write the sequences in the dataframe to a FASTA file.

    Args:
        df: DataFrame containing 'name' and 'sequence' columns.
        filename: Path to output FASTA file.
    """
    with open(filename, "w", encoding="utf-8") as f:
        for _, row in df.iterrows():
            f.write(f">{row['name']}\n")
            f.write(f"{row['sequence']}\n")


def to_opool(df: pd.DataFrame, name: str, filename: str) -> None:
    """Write the sequences in the dataframe to an opool file.

    Args:
        df: DataFrame containing sequences.
        name: Opool name to assign to all sequences.
        filename: Path to output opool (Excel) file.
    """
    df = df.copy()
    df["name"] = name
    df = df[["name", "sequence"]]
    df.to_xlsx(filename, index=False)


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


def trim(
    df: pd.DataFrame, start: int, end: int, extra_columns: list[str] = []
) -> pd.DataFrame:
    """
    Trims the 'sequence', 'structure', and 'data' columns of the DataFrame to the
    given start and end indices.

    Args:
        df (pd.DataFrame): A DataFrame with 'sequence', 'structure', and 'data'
                           columns, where 'data' contains lists of numbers.
        start (int): The start index for trimming.
        end (int): The end index for trimming.

    Returns:
        pd.DataFrame: A trimmed DataFrame with the 'sequence', 'structure', and
                      'data' columns adjusted to the specified indices.
    """

    def trim_column(column: pd.Series, start: int, end: int) -> pd.Series:
        if start == 0 and end == 0:
            return column
        if end == 0:
            return column.str[start:]
        elif start == 0:
            return column.str[:-end]
        else:
            return column.str[start:-end]

    df = df.copy()
    trim_columns = ["sequence"]
    if "structure" in df.columns:
        trim_columns.append("structure")
    trim_columns.extend(extra_columns)
    for col in trim_columns:
        if col in df.columns:
            if isinstance(df.iloc[0][col], list) or isinstance(
                df.iloc[0][col], np.ndarray
            ):
                df[col] = df[col].apply(
                    lambda x: x[start:-end] if end != 0 else x[start:]
                )
            else:
                df[col] = trim_column(df[col], start, end)

    return df


def trim_p5_and_p3(df: pd.DataFrame, extra_columns: list[str] = []) -> pd.DataFrame:
    """
    Trims the 5' and 3' ends of the data in the DataFrame.

    This function reads a CSV file containing p5 sequences, converts these
    sequences to RNA, checks for a common p5 sequence in the given DataFrame,
    and trims the DataFrame based on the length of this common p5 sequence and
    a fixed 3' end length.

    Args:
        df (pd.DataFrame): A DataFrame with a 'data' column containing
                           sequences as strings.

    Returns:
        pd.DataFrame: A trimmed DataFrame with the 5' and 3' ends trimmed.

    Raises:
        ValueError: If no common p5 sequence is found or the sequence is not
                    registered in the CSV file.
    """
    df_p5 = pd.read_csv(get_resources_path() / "p5_sequences.csv")
    df_p3 = pd.read_csv(get_resources_path() / "p3_sequences.csv")
    common_p5_seq = ""
    for p5_seq in df_p5["sequence"]:
        if has_5p_sequence(df, p5_seq):
            common_p5_seq = p5_seq
    for p3_seq in df_p3["sequence"]:
        if has_3p_sequence(df, p3_seq):
            common_p3_seq = p3_seq
    if len(common_p5_seq) == 0 or len(common_p3_seq) == 0:
        raise ValueError("No common p5 or p3 sequence found")
    return trim(df, len(common_p5_seq), len(common_p3_seq), extra_columns)


def transcribe(df: pd.DataFrame, ignore_missing_t7: bool = False) -> pd.DataFrame:
    """Transcribe DNA template sequences to RNA and remove T7 promoter.

    Converts DNA to RNA, removes the T7 promoter (first 20 bases), and folds
    the resulting RNA sequences.

    Args:
        df: DataFrame with DNA template sequences.
        ignore_missing_t7: If True, proceed even if sequences don't have T7 promoter.

    Returns:
        DataFrame with transcribed RNA sequences and folded structures.

    Raises:
        ValueError: If not all sequences have T7 promoter and ignore_missing_t7 is False.
    """
    if not has_t7_promoter(df) and not ignore_missing_t7:
        raise ValueError("not all sequences start with T7 promoter")
    if not ignore_missing_t7:
        df = trim(df, 20, 0)
    df = to_rna(df)
    df = fold(df)
    return df


def generate_mutated_sequences(
    template: str,
    num_mutations: int,
    num_sequences: int,
    p5_seq: Optional[str] = None,
    p3_seq: Optional[str] = None,
    ntype: str = "DNA",
) -> pd.DataFrame:
    """Generate mutated sequences from a template with optional constant 5' and 3' ends.

    Mutations are randomly distributed across the variable region (excluding
    constant 5' and 3' sequences if provided).

    Args:
        template: Template sequence to mutate.
        num_mutations: Number of mutations to introduce per sequence.
        num_sequences: Number of mutated sequences to generate.
        p5_seq: Optional constant 5' sequence (mutations only in middle region if provided).
        p3_seq: Optional constant 3' sequence (mutations only in middle region if provided).
        ntype: Nucleotide type: "DNA" or "RNA". Defaults to "DNA".

    Returns:
        DataFrame with mutated sequences in 'name' and 'sequence' columns.

    Raises:
        ValueError: If template is too short for specified constant sequences or mutations.
    """
    nucleotides = {"DNA": DEFAULT_DNA_NTS, "RNA": DEFAULT_RNA_NTS}
    nucs = nucleotides.get(ntype, DEFAULT_DNA_NTS)

    # Determine the variable region to mutate
    if p5_seq and p3_seq:
        p5_len = len(p5_seq)
        p3_len = len(p3_seq)
        if len(template) < p5_len + p3_len:
            raise ValueError(
                f"Template sequence ({len(template)} bp) is shorter than "
                f"combined 5' and 3' sequences ({p5_len + p3_len} bp)"
            )
        variable_region = template[p5_len:-p3_len] if p3_len > 0 else template[p5_len:]
    elif p5_seq:
        p5_len = len(p5_seq)
        if len(template) < p5_len:
            raise ValueError(
                f"Template sequence ({len(template)} bp) is shorter than "
                f"5' sequence ({p5_len} bp)"
            )
        variable_region = template[p5_len:]
    elif p3_seq:
        p3_len = len(p3_seq)
        if len(template) < p3_len:
            raise ValueError(
                f"Template sequence ({len(template)} bp) is shorter than "
                f"3' sequence ({p3_len} bp)"
            )
        variable_region = template[:-p3_len] if p3_len > 0 else template
    else:
        variable_region = template

    if len(variable_region) < num_mutations:
        raise ValueError(
            f"Variable region ({len(variable_region)} bp) is shorter than "
            f"number of mutations ({num_mutations})"
        )

    sequences = []
    names = []

    for i in range(num_sequences):
        # Create a mutable copy of the variable region
        mutated_var = list(variable_region)

        # Randomly select positions to mutate (without replacement)
        positions_to_mutate = random.sample(
            range(len(mutated_var)), min(num_mutations, len(mutated_var))
        )

        # Mutate each selected position
        for pos in positions_to_mutate:
            original_nuc = mutated_var[pos]
            # Choose a different nucleotide
            available_nucs = [n for n in nucs if n != original_nuc]
            mutated_var[pos] = random.choice(available_nucs)

        # Reconstruct the full sequence
        mutated_var_str = "".join(mutated_var)

        if p5_seq and p3_seq:
            full_sequence = p5_seq + mutated_var_str + p3_seq
        elif p5_seq:
            full_sequence = p5_seq + mutated_var_str
        elif p3_seq:
            full_sequence = mutated_var_str + p3_seq
        else:
            full_sequence = mutated_var_str

        sequences.append(full_sequence)
        names.append(f"mutated_seq_{i+1}")

    df = pd.DataFrame({"name": names, "sequence": sequences})
    return df


def generate_random_sequences(
    length: int,
    num_sequences: int,
    p5_seq: Optional[str] = None,
    p3_seq: Optional[str] = None,
    ntype: str = "DNA",
) -> pd.DataFrame:
    """Generate random sequences with optional constant 5' and 3' ends.

    The middle region between constant sequences (if provided) is randomly generated.

    Args:
        length: Total length of sequences (including constant 5' and 3' if provided).
        num_sequences: Number of random sequences to generate.
        p5_seq: Optional constant 5' sequence.
        p3_seq: Optional constant 3' sequence.
        ntype: Nucleotide type: "DNA" or "RNA". Defaults to "DNA".

    Returns:
        DataFrame with random sequences in 'name' and 'sequence' columns.

    Raises:
        ValueError: If length is not greater than the sum of constant sequence lengths.
    """
    nucleotides = {"DNA": DEFAULT_DNA_NTS, "RNA": DEFAULT_RNA_NTS}
    nucs = nucleotides.get(ntype, DEFAULT_DNA_NTS)

    # Calculate the length of the random middle region
    p5_len = len(p5_seq) if p5_seq else 0
    p3_len = len(p3_seq) if p3_seq else 0
    random_length = length - p5_len - p3_len

    if random_length <= 0:
        raise ValueError(
            f"Sequence length ({length}) must be greater than the sum of "
            f"5' ({p5_len}) and 3' ({p3_len}) constant sequences"
        )

    sequences = []
    names = []

    for i in range(num_sequences):
        # Generate random middle region
        random_seq = "".join(random.choice(nucs) for _ in range(random_length))

        # Reconstruct the full sequence
        if p5_seq and p3_seq:
            full_sequence = p5_seq + random_seq + p3_seq
        elif p5_seq:
            full_sequence = p5_seq + random_seq
        elif p3_seq:
            full_sequence = random_seq + p3_seq
        else:
            full_sequence = random_seq

        sequences.append(full_sequence)
        names.append(f"random_seq_{i+1}")

    df = pd.DataFrame({"name": names, "sequence": sequences})
    return df
