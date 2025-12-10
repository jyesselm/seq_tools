"""
I/O commands for seq_tools CLI.
"""

import click

from seq_tools import dataframe
from seq_tools.cli.common import get_input_dataframe
from seq_tools.logger import get_logger, setup_applevel_logger


@click.command(help="generate fasta file from csv")
@click.argument("data")
@click.option(
    "-o",
    "--output",
    help="output file (default: output.fasta)",
    default="output.fasta",
)
def to_fasta(data: str, output: str) -> None:
    """Generate FASTA file from CSV.

    Args:
        data: Can be a sequence or a file path.
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("to_fasta")
    df = get_input_dataframe(data)
    log.info(f"Writing FASTA output to: {output}")
    dataframe.to_fasta(df, output)
    log.info(f"Successfully wrote {len(df)} sequences to {output}")


@click.command(help="generate oligo pool file from csv")
@click.argument("data")
@click.option("-n", "--name", help="name of the opool file", default="opool")
@click.option(
    "-o",
    "--output",
    help="output file (default: output.xlsx)",
    default="output.xlsx",
)
def to_opool(data: str, name: str, output: str) -> None:
    """Generate oligo pool file from CSV.

    Args:
        data: Can be a sequence or a file path.
        name: Name of the opool file.
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("to_opool")
    df = get_input_dataframe(data)
    log.info(f"Writing opool output to: {output}")
    dataframe.to_opool(df, name, output)
    log.info(f"Successfully wrote {len(df)} sequences to {output}")
