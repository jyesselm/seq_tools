"""
module for working with dataframes that contain nucleotide sequences
"""

import os
import pandas as pd
import numpy as np
from typing import List, Optional
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import random

import editdistance
import vienna

from seq_tools import sequence, extinction_coeff
from seq_tools.structure import SequenceStructure
from seq_tools.structure import find as find_seq_struct
from seq_tools.config import T7_PROMOTER, DEFAULT_DNA_NTS, DEFAULT_RNA_NTS
from seq_tools.validation import validate_dataframe, ensure_name_column


def split(df: pd.DataFrame, n_chunks: int) -> List[pd.DataFrame]:
    """
    Splits a DataFrame into multiple chunks.

    Args:
        df (pd.DataFrame): The DataFrame to be split.
        n_chunks (int): The number of chunks to split the DataFrame into.

    Returns:
        List[pd.DataFrame]: A list of DataFrames, each representing a chunk of the original DataFrame.
    """
    chunk_size = len(df) // n_chunks
    chunks = [df.iloc[i * chunk_size : (i + 1) * chunk_size] for i in range(n_chunks)]

    # Handle the last chunk in case the division isn't perfect
    if len(df) % n_chunks != 0:
        chunks.append(df.iloc[n_chunks * chunk_size :])

    return chunks


def run_in_parallel(df: pd.DataFrame, func, threads: int) -> pd.DataFrame:
    """
    Runs a function in parallel on chunks of a DataFrame using multiple threads.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        func (callable): The function to apply to each chunk of the DataFrame.
        threads (int): The number of threads to use for parallel processing.

    Returns:
        pd.DataFrame: The combined results of applying the function to each chunk.
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
    """
    adds a 5' and 3' sequence to the sequences in the dataframe
    :param df: dataframe
    :param p5_seq: 5' sequence
    :param p3_seq: 3' sequence
    :return: None
    """
    df = df.copy()
    df["sequence"] = df["sequence"].apply(lambda x: p5_seq + x + p3_seq)
    if "structure" in df.columns:
        df = fold(df)
    return df


def calc_edit_distance(df: pd.DataFrame) -> float:
    """
    calculates the edit distance between each sequence in the dataframe
    :param df: dataframe
    :return: the edit distance
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


def _compute_pairwise_edit_distances(chunk):
    """
    Helper function to compute pairwise edit distances for a chunk of sequence pairs.
    Used for parallel processing.

    :param chunk: list of tuples (i, j, seq1, seq2)
    :return: dictionary mapping indices to their minimum edit distances
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
    """
    calculates the edit distance between each sequence in the dataframe using parallel processing

    :param df: dataframe
    :param n_workers: number of workers to use. If None, uses the number of CPU cores
    :param use_threads: if True, use threads instead of processes
    :return: the edit distance
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
    """
    determines the nucleotide type of the sequences in the dataframe
    :param df: dataframe
    :return: nucleotide type, RNA or DNA
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
    """
    folds each sequence in the dataframe
    :param df: dataframe
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
    """
    calculates the extinction coefficient for each sequence in the dataframe
    :param df: dataframe
    :param ntype: nucleotide type, RNA or DNA
    :param double_stranded: is double stranded?
    :return: None
    """

    def compute_w_struct(row) -> float:
        """
        computes the extinction coefficient for a sequence with a structure
        :param row: dataframe row
        :return: extinction coefficient
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
    """
    calculates the length of each sequence in the dataframe
    :param df: dataframe
    :return: None
    """
    df = df.copy()
    df["length"] = df["sequence"].apply(len)
    return df


def get_molecular_weight(
    df: pd.DataFrame, ntype: str, double_stranded: bool = False
) -> pd.DataFrame:
    """
    :param df: pandas data frame
    :param ntype: nucleotide type, RNA or DNA
    :param double_stranded: is double stranded?
    :return: None
    """
    df = df.copy()
    df["mw"] = df["sequence"].apply(
        lambda x: sequence.get_molecular_weight(x, ntype, double_stranded)
    )
    return df


def get_default_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds names to dataframe, if not already present
    :param df: dataframe
    :return: dataframe with names
    """
    if "name" in df.columns:
        raise ValueError("Dataframe already has names")
    # add `name` column to dataframe in the form of `seq_1`, `seq_2`, etc.
    df = df.copy()
    df["name"] = df.index.map(lambda x: "seq_" + str(x))
    return df


def get_reverse_complement(df: pd.DataFrame, ntype: str) -> pd.DataFrame:
    """
    reverse complements each sequence in the dataframe
    :param df: dataframe
    :param ntype: nucleotide type, RNA or DNA
    :return: stores reverse complement in dataframe rev_comp column
    """
    df = df.copy()
    df["rev_comp"] = df["sequence"].apply(
        lambda x: sequence.get_reverse_complement(x, ntype)
    )
    return df


def has_5p_sequence(df: pd.DataFrame, p5_seq: str) -> bool:
    """
    checks to see if p5_seq is present in the 5' end of the sequence
    :param df: dataframe
    :return: True if 5' sequence is present, False otherwise
    """
    return df["sequence"].str.startswith(p5_seq).all()


def has_3p_sequence(df: pd.DataFrame, p3_seq: str) -> bool:
    """
    checks to see if p5_seq is present in the 3' end of the sequence
    :param df: dataframe
    :return: True if 3' sequence is present, False otherwise
    """
    return df["sequence"].str.endswith(p3_seq).all()


def has_sequence(df: pd.DataFrame, seq: str) -> bool:
    """
    checks to see if seq is present in the sequence
    :param df: dataframe
    :return: True if sequence is present, False otherwise
    """
    return df["sequence"].str.contains(seq).all()


def has_t7_promoter(df: pd.DataFrame) -> bool:
    """
    Check if each sequence in the dataframe has a T7 promoter.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with sequences.

    Returns
    -------
    bool
        True if all sequences have T7 promoter, False otherwise.
    """
    has_t7 = df[df["sequence"].str.startswith(T7_PROMOTER)]
    if len(has_t7) != len(df):
        return False
    return True


def has_seq_struct(df: pd.DataFrame, seq_struct: SequenceStructure) -> bool:
    for _, row in df.iterrows():
        row_seq_struct = SequenceStructure(row["sequence"], row["structure"])
        if find_seq_struct(row_seq_struct, seq_struct) == 0:
            return False
    return True


def to_dna(df: pd.DataFrame) -> pd.DataFrame:
    """
    converts each sequence in dataframe to DNA
    :return: None
    """
    df = df.copy()
    df["sequence"] = df["sequence"].apply(sequence.to_dna)
    if "structure" in df.columns:
        df = df.drop(columns=["structure"])
    return df


def to_dna_template(df: pd.DataFrame) -> pd.DataFrame:
    """
    converts each sequence in dataframe to DNA
    :return: None
    """
    df = df.copy()
    df["sequence"] = df["sequence"].apply(sequence.to_dna_template)
    if "structure" in df.columns:
        df = df.drop(columns=["structure"])
    return df


def to_fasta(df: pd.DataFrame, filename: str) -> None:
    """
    writes the sequences in the dataframe to a fasta file
    :param df: dataframe
    :param filename: fasta file path
    :return: None
    """
    with open(filename, "w", encoding="utf-8") as f:
        for _, row in df.iterrows():
            f.write(f">{row['name']}\n")
            f.write(f"{row['sequence']}\n")


def to_opool(df: pd.DataFrame, name: str, filename: str) -> None:
    """
    writes the sequences in the dataframe to an opool file
    :param df: dataframe
    :param name: opool name
    :param filename: opool file path
    :return: None
    """
    df = df.copy()
    df["name"] = name
    df = df[["name", "sequence"]]
    df.to_xlsx(filename, index=False)


def to_rna(df: pd.DataFrame) -> pd.DataFrame:
    """
    converts each sequence in dataframe to DNA
    :return: None
    """
    df = df.copy()
    df["sequence"] = df["sequence"].apply(sequence.to_rna)
    return df


def trim(df: pd.DataFrame, p5_length: int, p3_length: int) -> pd.DataFrame:
    """
    takes a data frame and trims the sequences. If there is a structure
    it will also trim the structure
    :param df: dataframe
    :param p5_length: length to trim from 5'
    :param p3_length: length to trim from 3'
    :return: None
    """
    df = df.copy()
    # trim `sequence` column and `structure` column
    p3_length = -p3_length
    if p5_length == 0:
        p5_length = None
    if p3_length == 0:
        p3_length = None
    df["sequence"] = df["sequence"].str.slice(p5_length, p3_length)
    if "structure" in df.columns:
        df["structure"] = df["structure"].str.slice(p5_length, p3_length)
    return df


def transcribe(df: pd.DataFrame, ignore_missing_t7=False) -> pd.DataFrame:
    """
    transcribes each sequence in the dataframe (DNA -> RNA) removes t7 promoter
    :param df: dataframe with DNA template sequences
    :param ignore_missing_t7: ignore sequences that don't have a T7 promoter
    :return: dataframe with RNA sequences
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
    """
    Generate mutated sequences from a template sequence with optional constant 5' and 3' ends.

    Parameters
    ----------
    template : str
        Template sequence to mutate.
    num_mutations : int
        Number of mutations to introduce per sequence.
    num_sequences : int
        Number of mutated sequences to generate.
    p5_seq : str, optional
        Optional constant 5' sequence (if provided, mutations only in middle region).
    p3_seq : str, optional
        Optional constant 3' sequence (if provided, mutations only in middle region).
    ntype : str, optional
        Nucleotide type: "DNA" or "RNA" (default: "DNA").

    Returns
    -------
    pd.DataFrame
        DataFrame with mutated sequences in 'name' and 'sequence' columns.
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
    """
    Generate random sequences with optional constant 5' and 3' ends.

    Parameters
    ----------
    length : int
        Total length of sequences (including constant 5' and 3' if provided).
    num_sequences : int
        Number of random sequences to generate.
    p5_seq : str, optional
        Optional constant 5' sequence.
    p3_seq : str, optional
        Optional constant 3' sequence.
    ntype : str, optional
        Nucleotide type: "DNA" or "RNA" (default: "DNA").

    Returns
    -------
    pd.DataFrame
        DataFrame with random sequences in 'name' and 'sequence' columns.
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
