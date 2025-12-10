"""
Primer commands for seq_tools CLI.
"""

from typing import Optional

import click
import pandas as pd
import tabulate

from seq_tools import dataframe
from seq_tools.cli.common import get_input_dataframe, get_ntype, handle_output
from seq_tools.logger import get_logger, setup_applevel_logger
from seq_tools.utils import get_resources_path


def find_matching_sequences(
    target_seq: str, df_known: pd.DataFrame, check_start: bool = True
) -> list[pd.Series]:
    """Find known sequences that match the target sequence.

    Args:
        target_seq: Sequence to match against.
        df_known: DataFrame with known sequences.
        check_start: If True, check prefix matching; if False, check suffix.

    Returns:
        List of matching rows from df_known.
    """
    matching_known = []
    for _, row in df_known.iterrows():
        known_seq = row["sequence"]
        if check_start:
            if target_seq.startswith(known_seq) or known_seq.startswith(target_seq):
                matching_known.append(row)
        else:
            if target_seq.endswith(known_seq) or known_seq.endswith(target_seq):
                matching_known.append(row)
    return matching_known


def log_matching_sequences(matching_known: list[pd.Series], logger) -> None:
    """Log information about matching known sequences.

    Args:
        matching_known: List of matching sequences.
        logger: Logger instance.
    """
    if not matching_known:
        return

    logger.info("\nMatching known sequences:")
    for row in matching_known:
        name = row.get("name", "unknown")
        seq = row.get("sequence", "")
        code = row.get("code", "")
        structure = row.get("structure", "")
        logger.info(f"  - {name}: {seq}")
        if code:
            logger.info(f"    Code: {code}")
        if structure:
            logger.info(f"    Structure: {structure}")


@click.command(help="checks to see if p5 is present in all sequences")
@click.argument("data")
@click.option(
    "-p5",
    "--p5-seq",
    help="p5 sequence (if not provided, finds longest common prefix)",
    default=None,
)
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
def has_p5(data: str, p5_seq: Optional[str], ntype: Optional[str]) -> None:
    """Check if a sequence has a p5 sequence.

    If no p5 sequence is provided, automatically finds the longest
    common prefix from all sequences.

    Args:
        data: Can be a sequence or a file path.
        p5_seq: P5 sequence (optional).
        ntype: Type of nucleic acid.
    """
    setup_applevel_logger()
    log = get_logger("has_p5")
    df = get_input_dataframe(data)

    if ntype is not None:
        get_ntype(df, ntype)

    df_p5 = pd.read_csv(get_resources_path() / "p5_sequences.csv")

    if p5_seq is None:
        common_prefix = dataframe.find_longest_common_prefix(df)
        if common_prefix:
            log.info(
                f"Found longest common 5' sequence: {common_prefix} "
                f"(length: {len(common_prefix)})"
            )
            log.info("This sequence is present at the 5' end of all sequences")
            matching_known = find_matching_sequences(
                common_prefix, df_p5, check_start=True
            )
            if matching_known:
                log_matching_sequences(matching_known, log)
            else:
                log.info("No matching known p5 sequences found in resource file")
        else:
            log.info("No common 5' sequence found in all sequences")
    else:
        has_p5_seq = dataframe.has_5p_sequence(df, p5_seq)
        if has_p5_seq:
            log.info(f"p5 sequence '{p5_seq}' is present in all sequences")
            matching_known = find_matching_sequences(p5_seq, df_p5, check_start=True)
            if matching_known:
                log_matching_sequences(matching_known, log)
        else:
            log.info(f"p5 sequence '{p5_seq}' is not present in all sequences")


@click.command(help="checks to see if p3 is present in all sequences")
@click.argument("data")
@click.option(
    "-p3",
    "--p3-seq",
    help="p3 sequence (if not provided, finds longest common suffix)",
    default=None,
)
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
def has_p3(data: str, p3_seq: Optional[str], ntype: Optional[str]) -> None:
    """Check if a sequence has a p3 sequence.

    If no p3 sequence is provided, automatically finds the longest
    common suffix from all sequences.

    Args:
        data: Can be a sequence or a file path.
        p3_seq: P3 sequence (optional).
        ntype: Type of nucleic acid.
    """
    setup_applevel_logger()
    log = get_logger("has_p3")
    df = get_input_dataframe(data)

    if ntype is not None:
        get_ntype(df, ntype)

    df_p3 = pd.read_csv(get_resources_path() / "p3_sequences.csv")

    if p3_seq is None:
        common_suffix = dataframe.find_longest_common_suffix(df)
        if common_suffix:
            log.info(
                f"Found longest common 3' sequence: {common_suffix} "
                f"(length: {len(common_suffix)})"
            )
            log.info("This sequence is present at the 3' end of all sequences")
            matching_known = find_matching_sequences(
                common_suffix, df_p3, check_start=False
            )
            if matching_known:
                log_matching_sequences(matching_known, log)
            else:
                log.info("No matching known p3 sequences found in resource file")
        else:
            log.info("No common 3' sequence found in all sequences")
    else:
        has_p3_seq = dataframe.has_3p_sequence(df, p3_seq)
        if has_p3_seq:
            log.info(f"p3 sequence '{p3_seq}' is present in all sequences")
            matching_known = find_matching_sequences(p3_seq, df_p3, check_start=False)
            if matching_known:
                log_matching_sequences(matching_known, log)
        else:
            log.info(f"p3 sequence '{p3_seq}' is not present in all sequences")


def resolve_p5_sequence(
    p5_name: Optional[str],
    p5_seq: Optional[str],
    p5_ss: Optional[str],
    df_p5: pd.DataFrame,
) -> tuple[str, str]:
    """Resolve p5 sequence and structure from name or custom values.

    Args:
        p5_name: Name of p5 sequence from resources.
        p5_seq: Custom p5 sequence.
        p5_ss: Custom p5 structure.
        df_p5: DataFrame with p5 sequences from resources.

    Returns:
        Tuple of (sequence, structure).
    """
    if p5_name:
        p5_row = df_p5[df_p5["name"] == p5_name]
        if p5_row.empty:
            raise click.ClickException(
                f"P5 sequence '{p5_name}' not found in resources. "
                "Use 'seq-tools list-common-seqs' to see available sequences."
            )
        return p5_row.iloc[0]["sequence"], p5_row.iloc[0]["structure"]

    if p5_seq:
        if not p5_ss:
            raise click.ClickException(
                "When using --p5-seq, you must also provide --p5-ss"
            )
        return p5_seq, p5_ss

    return "", ""


def resolve_p3_sequence(
    p3_name: Optional[str],
    p3_seq: Optional[str],
    p3_ss: Optional[str],
    df_p3: pd.DataFrame,
) -> tuple[str, str]:
    """Resolve p3 sequence and structure from name or custom values.

    Args:
        p3_name: Name of p3 sequence from resources.
        p3_seq: Custom p3 sequence.
        p3_ss: Custom p3 structure.
        df_p3: DataFrame with p3 sequences from resources.

    Returns:
        Tuple of (sequence, structure).
    """
    if p3_name:
        p3_row = df_p3[df_p3["name"] == p3_name]
        if p3_row.empty:
            raise click.ClickException(
                f"P3 sequence '{p3_name}' not found in resources. "
                "Use 'seq-tools list-common-seqs' to see available sequences."
            )
        return p3_row.iloc[0]["sequence"], p3_row.iloc[0]["structure"]

    if p3_seq:
        if not p3_ss:
            raise click.ClickException(
                "When using --p3-seq, you must also provide --p3-ss"
            )
        return p3_seq, p3_ss

    return "", ""


@click.command(
    name="add-common-seqs",
    help="add common 5' and/or 3' sequences and validate structure",
)
@click.argument("data")
@click.option("-p5", "--p5-name", help="name of 5' common sequence from resources")
@click.option("-p3", "--p3-name", help="name of 3' common sequence from resources")
@click.option("--p5-seq", help="custom 5' sequence")
@click.option("--p5-ss", help="custom 5' secondary structure (dot-bracket)")
@click.option("--p3-seq", help="custom 3' sequence")
@click.option("--p3-ss", help="custom 3' secondary structure (dot-bracket)")
@click.option(
    "-p",
    "--parallel",
    is_flag=True,
    help="use parallel processing for structure validation",
)
@click.option(
    "-w",
    "--workers",
    type=int,
    help="number of workers for parallel processing",
    default=4,
)
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def add_common_seqs(
    data: str,
    p5_name: Optional[str],
    p3_name: Optional[str],
    p5_seq: Optional[str],
    p5_ss: Optional[str],
    p3_seq: Optional[str],
    p3_ss: Optional[str],
    parallel: bool,
    workers: int,
    output: str,
) -> None:
    """Add common 5' and/or 3' sequences and validate structure.

    You can either:
    1. Specify custom sequences and structures with --p5-seq/--p5-ss
       and/or --p3-seq/--p3-ss
    2. Use names from the resources folder with -p5/-p3

    Args:
        data: Input CSV file with sequences and structures.
        p5_name: Name of 5' common sequence from resources.
        p3_name: Name of 3' common sequence from resources.
        p5_seq: Custom 5' sequence.
        p5_ss: Custom 5' structure.
        p3_seq: Custom 3' sequence.
        p3_ss: Custom 3' structure.
        parallel: Use parallel processing.
        workers: Number of workers.
        output: Output file.
    """
    setup_applevel_logger()
    log = get_logger("add_common_seqs")
    df = get_input_dataframe(data)

    if "structure" not in df.columns:
        raise click.ClickException("Input data must contain a 'structure' column")

    df_p5 = pd.read_csv(get_resources_path() / "p5_sequences.csv")
    df_p3 = pd.read_csv(get_resources_path() / "p3_sequences.csv")

    p5_sequence, p5_structure = resolve_p5_sequence(p5_name, p5_seq, p5_ss, df_p5)
    p3_sequence, p3_structure = resolve_p3_sequence(p3_name, p3_seq, p3_ss, df_p3)

    if p5_sequence:
        log.info(
            f"Using p5 '{p5_name or 'custom'}': {p5_sequence} "
            f"(structure: {p5_structure})"
        )
    if p3_sequence:
        log.info(
            f"Using p3 '{p3_name or 'custom'}': {p3_sequence} "
            f"(structure: {p3_structure})"
        )

    if not p5_sequence and not p3_sequence:
        raise click.ClickException(
            "You must specify either p5 or p3 sequences (or both)"
        )

    log.info(f"Processing {len(df)} sequences...")
    result_df = dataframe.add_common_seqs(
        df,
        p5_sequence,
        p5_structure,
        p3_sequence,
        p3_structure,
        parallel=parallel,
        workers=workers,
    )

    matches = result_df["structure_match"].sum()
    total = len(result_df)
    log.info(
        f"Structure validation: {matches}/{total} sequences match expected structure"
    )

    if matches < total:
        log.warning(f"{total - matches} sequences did not match expected structure")

    handle_output(result_df, output, show_all=True)


@click.command(
    name="list-common-seqs", help="list available 5' and 3' common sequences"
)
def list_common_seqs() -> None:
    """List all available 5' and 3' common sequences from resource files."""
    setup_applevel_logger()
    log = get_logger("list_common_seqs")

    p5_path = get_resources_path() / "p5_sequences.csv"
    if p5_path.exists():
        df_p5 = pd.read_csv(p5_path)
        log.info("\nAvailable 5' sequences:")
        log.info(
            "\n"
            + tabulate.tabulate(
                df_p5, headers="keys", tablefmt="simple", showindex=False
            )
        )
    else:
        log.warning(f"p5_sequences.csv not found at {p5_path}")

    p3_path = get_resources_path() / "p3_sequences.csv"
    if p3_path.exists():
        df_p3 = pd.read_csv(p3_path)
        log.info("\n\nAvailable 3' sequences:")
        log.info(
            "\n"
            + tabulate.tabulate(
                df_p3, headers="keys", tablefmt="simple", showindex=False
            )
        )
    else:
        log.warning(f"p3_sequences.csv not found at {p3_path}")


@click.command(
    name="remove-common-seqs",
    help="identify and remove common 5' and 3' sequences",
)
@click.argument("data")
@click.option(
    "-o",
    "--output",
    help="output file (default: output.csv)",
    default="output.csv",
)
def remove_common_seqs(data: str, output: str) -> None:
    """Identify and remove common 5' and 3' sequences from sequences.

    This command checks all sequences in the input CSV against known
    5' and 3' sequences from resource files. It attempts to match both
    sequence and structure patterns, and provides warnings if only
    sequence or only structure matches.

    Args:
        data: Input CSV file with sequences and optionally structures.
        output: Output CSV file.
    """
    setup_applevel_logger()
    log = get_logger("remove_common_seqs")
    df = get_input_dataframe(data)

    try:
        df, result = dataframe.remove_common_seqs(df)

        removed = result.get("removed", {})
        if removed.get("p5_sequence"):
            log.info(
                f"Removed 5' sequence: {removed['p5_sequence']} "
                f"(length: {removed['p5_length']})"
            )
        if removed.get("p3_sequence"):
            log.info(
                f"Removed 3' sequence: {removed['p3_sequence']} "
                f"(length: {removed['p3_length']})"
            )

        if result.get("warnings"):
            for warning in result["warnings"]:
                log.warning(warning)

        log.info("Successfully identified and removed common 5' and/or 3' sequences")

        if len(df) > 0:
            log.info(f"\nFirst sequence (showing 1 of {len(df)}):")
            log.info(
                "\n"
                + tabulate.tabulate(
                    df.iloc[0:1],
                    headers="keys",
                    tablefmt="simple",
                    showindex=False,
                )
            )
        df.to_csv(output, index=False)
        log.info(f"\nSuccessfully wrote {len(df)} sequences to {output}")
    except ValueError as e:
        log.error(f"Error: {e}")
        raise click.ClickException(str(e)) from e
