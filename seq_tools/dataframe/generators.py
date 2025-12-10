"""Sequence generation functions for dataframes."""

import random
from typing import Optional

import pandas as pd

from seq_tools.config import DEFAULT_DNA_NTS, DEFAULT_RNA_NTS


def _calculate_variable_region(
    template: str,
    p5_seq: Optional[str],
    p3_seq: Optional[str],
) -> tuple[str, int, int]:
    """Calculate the variable region that can be mutated.

    Args:
        template: Template sequence.
        p5_seq: Optional 5' constant sequence.
        p3_seq: Optional 3' constant sequence.

    Returns:
        Tuple of (variable_region, p5_len, p3_len).

    Raises:
        ValueError: If template is too short for specified constant sequences.
    """
    p5_len = len(p5_seq) if p5_seq else 0
    p3_len = len(p3_seq) if p3_seq else 0

    if p5_seq and p3_seq:
        if len(template) < p5_len + p3_len:
            raise ValueError(
                f"Template sequence ({len(template)} bp) is shorter than "
                f"combined 5' and 3' sequences ({p5_len + p3_len} bp)"
            )
        if p3_len > 0:
            variable_region = template[p5_len:-p3_len]
        else:
            variable_region = template[p5_len:]
    elif p5_seq:
        if len(template) < p5_len:
            raise ValueError(
                f"Template sequence ({len(template)} bp) is shorter than "
                f"5' sequence ({p5_len} bp)"
            )
        variable_region = template[p5_len:]
    elif p3_seq:
        if len(template) < p3_len:
            raise ValueError(
                f"Template sequence ({len(template)} bp) is shorter than "
                f"3' sequence ({p3_len} bp)"
            )
        if p3_len > 0:
            variable_region = template[:-p3_len]
        else:
            variable_region = template
    else:
        variable_region = template

    return variable_region, p5_len, p3_len


def _mutate_sequence(variable_region: str, num_mutations: int, nucleotides: str) -> str:
    """Mutate a variable region sequence.

    Args:
        variable_region: The region to mutate.
        num_mutations: Number of mutations to introduce.
        nucleotides: Available nucleotides for mutation.

    Returns:
        Mutated variable region sequence.
    """
    mutated_var = list(variable_region)

    # Randomly select positions to mutate (without replacement)
    positions_to_mutate = random.sample(
        range(len(mutated_var)), min(num_mutations, len(mutated_var))
    )

    # Mutate each selected position
    for pos in positions_to_mutate:
        original_nuc = mutated_var[pos]
        # Choose a different nucleotide
        available_nucs = [n for n in nucleotides if n != original_nuc]
        mutated_var[pos] = random.choice(available_nucs)

    return "".join(mutated_var)


def _reconstruct_full_sequence(
    mutated_var: str,
    p5_seq: Optional[str],
    p3_seq: Optional[str],
) -> str:
    """Reconstruct full sequence from mutated variable region.

    Args:
        mutated_var: Mutated variable region.
        p5_seq: Optional 5' constant sequence.
        p3_seq: Optional 3' constant sequence.

    Returns:
        Full reconstructed sequence.
    """
    if p5_seq and p3_seq:
        return p5_seq + mutated_var + p3_seq
    elif p5_seq:
        return p5_seq + mutated_var
    elif p3_seq:
        return mutated_var + p3_seq
    else:
        return mutated_var


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
    variable_region, _, _ = _calculate_variable_region(template, p5_seq, p3_seq)

    if len(variable_region) < num_mutations:
        raise ValueError(
            f"Variable region ({len(variable_region)} bp) is shorter than "
            f"number of mutations ({num_mutations})"
        )

    sequences = []
    names = []

    for i in range(num_sequences):
        # Mutate the variable region
        mutated_var = _mutate_sequence(variable_region, num_mutations, nucs)

        # Reconstruct the full sequence
        full_sequence = _reconstruct_full_sequence(mutated_var, p5_seq, p3_seq)

        sequences.append(full_sequence)
        names.append(f"mutated_seq_{i + 1}")

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
        full_sequence = _reconstruct_full_sequence(random_seq, p5_seq, p3_seq)

        sequences.append(full_sequence)
        names.append(f"random_seq_{i + 1}")

    df = pd.DataFrame({"name": names, "sequence": sequences})
    return df
