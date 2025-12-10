"""
Generator commands for seq_tools CLI.
"""

import os
from typing import Optional

import click
import pandas as pd

from seq_tools import dataframe
from seq_tools.cli.common import (
    handle_output,
    validate_dataframe,
)
from seq_tools.logger import get_logger, setup_applevel_logger


def get_template_sequence(template: str) -> str:
    """Get template sequence from file or string.

    Args:
        template: Template sequence or file path.

    Returns:
        Template sequence as string.
    """
    log = get_logger("get_template_sequence")

    if os.path.isfile(template):
        log.info(f"reading template from file {template}")
        df_temp = pd.read_csv(template)
        validate_dataframe(df_temp)
        if len(df_temp) != 1:
            raise ValueError(
                f"Template file must contain exactly one sequence, found {len(df_temp)}"
            )
        return df_temp.iloc[0]["sequence"]

    log.info(f"using template sequence: {template}")
    return template


@click.command(
    help="generate mutated sequences from a template with optional constant ends"
)
@click.argument("template")
@click.option(
    "-n",
    "--num-sequences",
    type=int,
    default=10,
    help="number of mutated sequences to generate",
)
@click.option(
    "-m",
    "--num-mutations",
    type=int,
    required=True,
    help="number of mutations per sequence",
)
@click.option(
    "-p5",
    "--p5-seq",
    default=None,
    help="constant 5' sequence (mutations only in middle region)",
)
@click.option(
    "-p3",
    "--p3-seq",
    default=None,
    help="constant 3' sequence (mutations only in middle region)",
)
@click.option(
    "-nt",
    "--ntype",
    default="DNA",
    type=click.Choice(["DNA", "RNA"]),
    help="type of nucleic acid",
)
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def mutate(
    template: str,
    num_sequences: int,
    num_mutations: int,
    p5_seq: Optional[str],
    p3_seq: Optional[str],
    ntype: str,
    output: str,
) -> None:
    """Generate mutated sequences from a template sequence.

    Optional constant 5' and 3' ends can be specified to restrict
    mutations to the middle region only.

    Args:
        template: Template sequence to mutate (can be a sequence or file).
        num_sequences: Number of mutated sequences to generate.
        num_mutations: Number of mutations to introduce per sequence.
        p5_seq: Optional constant 5' sequence.
        p3_seq: Optional constant 3' sequence.
        ntype: Type of nucleic acid (DNA or RNA).
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("mutate")

    template_seq = get_template_sequence(template)

    log.info(
        f"generating {num_sequences} sequences with {num_mutations} mutations each"
    )
    if p5_seq:
        log.info(f"using constant 5' sequence: {p5_seq}")
    if p3_seq:
        log.info(f"using constant 3' sequence: {p3_seq}")

    df = dataframe.generate_mutated_sequences(
        template=template_seq,
        num_mutations=num_mutations,
        num_sequences=num_sequences,
        p5_seq=p5_seq,
        p3_seq=p3_seq,
        ntype=ntype,
    )
    handle_output(df, output)


@click.command(help="generate random sequences with optional constant 5' and 3' ends")
@click.option(
    "-l",
    "--length",
    type=int,
    required=True,
    help="total length of sequences (including constant 5' and 3')",
)
@click.option(
    "-n",
    "--num-sequences",
    type=int,
    default=10,
    help="number of random sequences to generate",
)
@click.option(
    "-p5",
    "--p5-seq",
    default=None,
    help="constant 5' sequence",
)
@click.option(
    "-p3",
    "--p3-seq",
    default=None,
    help="constant 3' sequence",
)
@click.option(
    "-nt",
    "--ntype",
    default="DNA",
    type=click.Choice(["DNA", "RNA"]),
    help="type of nucleic acid",
)
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def random(
    length: int,
    num_sequences: int,
    p5_seq: Optional[str],
    p3_seq: Optional[str],
    ntype: str,
    output: str,
) -> None:
    """Generate random sequences with optional constant 5' and 3' ends.

    Args:
        length: Total length of sequences (including constant 5' and 3').
        num_sequences: Number of random sequences to generate.
        p5_seq: Optional constant 5' sequence.
        p3_seq: Optional constant 3' sequence.
        ntype: Type of nucleic acid (DNA or RNA).
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("random")

    log.info(f"generating {num_sequences} random sequences of length {length}")
    if p5_seq:
        log.info(f"using constant 5' sequence: {p5_seq}")
    if p3_seq:
        log.info(f"using constant 3' sequence: {p3_seq}")

    df = dataframe.generate_random_sequences(
        length=length,
        num_sequences=num_sequences,
        p5_seq=p5_seq,
        p3_seq=p3_seq,
        ntype=ntype,
    )
    handle_output(df, output)
