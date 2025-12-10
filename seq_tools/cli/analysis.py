"""
Analysis commands for seq_tools CLI.
"""

from typing import Optional

import click
import pandas as pd

from seq_tools import dataframe
from seq_tools.cli.common import get_input_dataframe, get_ntype, handle_output
from seq_tools.logger import get_logger, setup_applevel_logger


@click.command(help="calculate the edit distance of a library")
@click.option(
    "-p",
    "--parallel",
    is_flag=True,
    help="calculate the edit distance in parallel",
)
@click.option("-w", "--workers", type=int, help="number of workers to use", default=1)
@click.option(
    "--use-threads",
    is_flag=True,
    help="use threads instead of processes",
    default=False,
)
@click.argument("data", type=click.Path(exists=True))
def edit_distance(data: str, parallel: bool, workers: int, use_threads: bool) -> None:
    """Calculate the edit distance of a library.

    Args:
        data: File path to CSV with sequences.
        parallel: Whether to use parallel processing.
        workers: Number of workers for parallel processing.
        use_threads: Whether to use threads instead of processes.
    """
    setup_applevel_logger()
    log = get_logger("edit_distance")
    df = pd.read_csv(data)

    if parallel:
        log.info(f"calculating edit distance in parallel with {workers} workers")
        if use_threads:
            log.info("using threads instead of processes")
        score = dataframe.calc_edit_distance_parallel(df, workers, use_threads)
    else:
        score = dataframe.calc_edit_distance(df)

    log.info(f"edit distance: {score}")


@click.command(help="calculate the extinction coefficient for each sequence")
@click.argument("data")
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
@click.option("-ds", "--double-stranded", is_flag=True)
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def ec(data: str, ntype: Optional[str], double_stranded: bool, output: str) -> None:
    """Calculate the extinction coefficient for each sequence.

    Args:
        data: Can be a sequence or a file path.
        ntype: Type of nucleic acid (RNA, DNA, or None).
        double_stranded: Whether sequence is double stranded.
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("extinction_coeff")
    df = get_input_dataframe(data)
    ntype = get_ntype(df, ntype)
    df = dataframe.get_extinction_coeff(df, ntype, double_stranded)

    if len(df) > 0 and "extinction_coeff" in df.columns:
        first_ec = df.iloc[0]["extinction_coeff"]
        strand_type = "double-stranded" if double_stranded else "single-stranded"
        log.info(f"Calculated extinction coefficients for {ntype} ({strand_type})")
        if len(df) == 1:
            log.info(f"Extinction coefficient: {first_ec:.0f} M⁻¹cm⁻¹")
        else:
            log.info(f"First sequence extinction coefficient: {first_ec:.0f} M⁻¹cm⁻¹")
            log.info(
                f"Average extinction coefficient: "
                f"{df['extinction_coeff'].mean():.0f} M⁻¹cm⁻¹"
            )

    handle_output(df, output)


@click.command(help="calculate the molecular weight for each sequence")
@click.argument("data")
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
@click.option("-ds", "--double-stranded", is_flag=True)
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def mw(data: str, ntype: Optional[str], double_stranded: bool, output: str) -> None:
    """Calculate the molecular weight for each sequence.

    Args:
        data: Can be a sequence or a file path.
        ntype: Type of nucleic acid (RNA, DNA, or None).
        double_stranded: Whether sequence is double stranded.
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("molecular_weight")
    df = get_input_dataframe(data)
    ntype = get_ntype(df, ntype)
    df = dataframe.get_molecular_weight(df, ntype, double_stranded)

    if len(df) > 0 and "mw" in df.columns:
        first_mw = df.iloc[0]["mw"]
        strand_type = "double-stranded" if double_stranded else "single-stranded"
        log.info(f"Calculated molecular weights for {ntype} ({strand_type})")
        if len(df) == 1:
            log.info(f"Molecular weight: {first_mw:.2f} Da")
        else:
            log.info(f"First sequence molecular weight: {first_mw:.2f} Da")
            log.info(f"Average molecular weight: {df['mw'].mean():.2f} Da")

    handle_output(df, output)


@click.command(help="calculate reverse complement for each sequence")
@click.argument("data")
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def rc(data: str, ntype: Optional[str], output: str) -> None:
    """Calculate the reverse complement for each sequence.

    Args:
        data: Can be a sequence or a file path.
        ntype: Type of nucleic acid (RNA, DNA, or None).
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("rc")
    df = get_input_dataframe(data)
    ntype = get_ntype(df, ntype)
    original_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
    df = dataframe.get_reverse_complement(df, ntype)

    if original_seq and len(df) > 0:
        rev_comp = df.iloc[0]["rev_comp"] if "rev_comp" in df.columns else ""
        log.info(f"Calculated reverse complement for {ntype}")
        if len(df) == 1:
            log.info(f"Example: {original_seq[:20]}... -> {rev_comp[:20]}...")

    handle_output(df, output)
