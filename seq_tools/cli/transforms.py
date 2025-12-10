"""
Transform commands for seq_tools CLI.
"""

import click

from seq_tools import dataframe, sequence
from seq_tools.cli.common import get_input_dataframe, handle_output
from seq_tools.logger import get_logger, setup_applevel_logger


@click.command(help="add a sequence to 5' and/or 3'")
@click.argument("data")
@click.option("-p5", "--p5-seq", default="")
@click.option("-p3", "--p3-seq", default="")
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def add(data: str, p5_seq: str, p3_seq: str, output: str) -> None:
    """Add a sequence to a dataframe at 5' and/or 3' ends.

    Args:
        data: Can be a sequence or a file path.
        p5_seq: Sequence to add to 5' end.
        p3_seq: Sequence to add to 3' end.
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("add")
    df = get_input_dataframe(data)

    has_structure = "structure" in df.columns

    if p5_seq:
        log.info(f"Adding 5' sequence: {p5_seq} (length: {len(p5_seq)})")
    if p3_seq:
        log.info(f"Adding 3' sequence: {p3_seq} (length: {len(p3_seq)})")

    if has_structure:
        log.info(
            "Structure column detected - sequences will be refolded "
            "after adding primers"
        )

    df = dataframe.add(df, p5_seq, p3_seq)
    handle_output(df, output)


@click.command(help="fold rna sequences")
@click.argument("data")
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def fold(data: str, output: str) -> None:
    """Fold RNA sequences using ViennaRNA.

    Args:
        data: Can be a sequence or a file path.
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("fold")
    df = get_input_dataframe(data)
    df = dataframe.fold(df)

    if len(df) > 0 and "structure" in df.columns:
        first_row = df.iloc[0]
        log.info("Folded sequences using ViennaRNA")
        if len(df) == 1:
            log.info(f"Structure: {first_row.get('structure', 'N/A')}")
            log.info(f"MFE: {first_row.get('mfe', 'N/A')} kcal/mol")
            log.info(f"Ensemble defect: {first_row.get('ens_defect', 'N/A')}")

    handle_output(df, output)


@click.command(help="convert rna sequence(s) to dna")
@click.argument("data")
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def to_dna(data: str, output: str) -> None:
    """Convert RNA sequence to DNA (U -> T).

    Args:
        data: Can be a sequence or a file path.
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("to_dna")
    df = get_input_dataframe(data)
    original_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
    df = df[["name", "sequence"]]
    df = dataframe.to_dna(df)

    if original_seq:
        converted_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
        log.info("Converted RNA to DNA (U -> T)")
        if len(df) == 1:
            log.info(f"Example: {original_seq[:20]}... -> {converted_seq[:20]}...")

    handle_output(df, output)


@click.command(help="convert rna sequence(s) to dna template, includes T7 promoter")
@click.argument("data")
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def to_dna_template(data: str, output: str) -> None:
    """Convert RNA sequence to DNA template with T7 promoter.

    Args:
        data: Can be a sequence or a file path.
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("to_dna_template")
    df = get_input_dataframe(data)
    original_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
    df = df[["name", "sequence"]]
    df = dataframe.to_dna_template(df)

    if original_seq:
        converted_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
        log.info("Converted RNA to DNA template with T7 promoter")
        if len(df) == 1:
            log.info("Added T7 promoter (20 bases) and converted U -> T")
            log.info(f"Example: {original_seq[:20]}... -> {converted_seq[:40]}...")

    handle_output(df, output)


@click.command(help="convert dna sequence(s) to rna")
@click.argument("data")
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def to_rna(data: str, output: str) -> None:
    """Convert DNA sequence to RNA (T -> U).

    Args:
        data: Can be a sequence or a file path.
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("to_rna")
    df = get_input_dataframe(data)
    original_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
    df = df[["name", "sequence"]]
    df["sequence"] = df["sequence"].apply(sequence.to_rna)

    if original_seq:
        converted_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
        log.info("Converted DNA to RNA (T -> U)")
        if len(df) == 1:
            log.info(f"Example: {original_seq[:20]}... -> {converted_seq[:20]}...")

    handle_output(df, output)


@click.command(help="trim 5'/3' ends of sequences")
@click.argument("data")
@click.option("-p5", "--p5-cut", default=0)
@click.option("-p3", "--p3-cut", default=0)
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def trim(data: str, p5_cut: int, p3_cut: int, output: str) -> None:
    """Trim 5' and/or 3' ends of sequences.

    Args:
        data: Can be a sequence or a file path.
        p5_cut: Number of bases to trim from 5' end.
        p3_cut: Number of bases to trim from 3' end.
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("trim")
    df = get_input_dataframe(data)

    if p5_cut > 0:
        log.info(f"Trimming {p5_cut} bases from 5' end")
    if p3_cut > 0:
        log.info(f"Trimming {p3_cut} bases from 3' end")

    df = dataframe.trim(df, p5_cut, p3_cut)
    handle_output(df, output)


@click.command(help="convert dna sequence(s) to rna")
@click.argument("data")
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def transcribe(data: str, output: str) -> None:
    """Transcribe DNA sequence to RNA (removes T7 promoter if present).

    Args:
        data: Can be a sequence or a file path.
        output: Output file path.
    """
    setup_applevel_logger()
    log = get_logger("transcribe")
    df = get_input_dataframe(data)
    original_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
    df = df[["name", "sequence"]]
    df = dataframe.transcribe(df)

    if original_seq:
        transcribed_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
        log.info("Transcribed DNA to RNA (removed T7 promoter if present, T -> U)")
        if len(df) == 1:
            log.info(f"Example: {original_seq[:30]}... -> {transcribed_seq[:20]}...")

    handle_output(df, output)
