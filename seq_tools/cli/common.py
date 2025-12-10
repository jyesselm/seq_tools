"""
Common utilities for CLI commands.
"""

import os
from typing import Optional

import pandas as pd
import tabulate

from seq_tools.logger import get_logger
from seq_tools.validation import ensure_name_column
from seq_tools.validation import validate_dataframe as validate_df


def validate_dataframe(df: pd.DataFrame) -> None:
    """Validate a dataframe to have a column named `sequence` and `name`.

    Args:
        df: Dataframe with sequences.
    """
    validate_df(df, require_name=False)
    ensure_name_column(df)


def get_input_dataframe(data: str) -> pd.DataFrame:
    """Return a dataframe from a sequence or a file.

    Args:
        data: Can be a sequence or a file path.

    Returns:
        DataFrame with sequences.
    """
    log = get_logger("get_input_dataframe")
    if os.path.isfile(data):
        log.info(f"reading file {data}")
        df = pd.read_csv(data)
        log.info(f"csv file contains {len(df)} sequences")
    else:
        log.info(f"reading sequence {data}")
        data_df = [["seq", data]]
        df = pd.DataFrame(data_df, columns=["name", "sequence"])
    validate_dataframe(df)
    return df


def get_ntype(df: pd.DataFrame, ntype: Optional[str]) -> str:
    """Handle the ntype parameter and determine nucleotide type.

    Args:
        df: Dataframe with sequences.
        ntype: Nucleotide type (DNA, RNA, or None).

    Returns:
        Nucleotide type as string.
    """
    from seq_tools import dataframe

    log = get_logger("handle_ntype")
    df_ntype = dataframe.determine_ntype(df)
    log.info(f"determining nucleic acid type: {df_ntype}")

    if ntype == "DNA":
        log.info("forcing sequences to be DNA")
        dataframe.to_dna(df)
        return ntype
    if ntype == "RNA":
        log.info("forcing sequences to be RNA")
        dataframe.to_rna(df)
        return ntype
    return df_ntype


def format_series_output(series: pd.Series) -> str:
    """Format a pandas Series for clean output display matching pandas format.

    Pandas aligns all lines to the same total length:
    total_length = max_key_length + max_value_length + 4

    Args:
        series: Pandas Series to format.

    Returns:
        Formatted string representation.
    """
    # Use pandas' native string formatting which handles floats nicely
    output = str(series)
    # Remove the dtype line at the end
    lines = output.splitlines()[:-1]
    return "\n".join(lines)


def handle_output(df: pd.DataFrame, output: str, show_all: bool = False) -> None:
    """Handle the output of the dataframe.

    Args:
        df: Dataframe with sequences.
        output: Output file path.
        show_all: If True, show all sequences; if False, show only first.
    """
    log = get_logger("handle_output")
    if len(df) == 1:
        log.info(f"output->\n{format_series_output(df.iloc[0])}")
        log.info(f"Writing output to: {output}")
        df.to_csv(output, index=False)
    else:
        log.info(f"Writing output CSV to: {output}")
        if show_all:
            if len(df) > 100:
                log.info(
                    "\n"
                    + tabulate.tabulate(df[0:100], headers="keys", tablefmt="simple")
                )
                log.info(f"... (showing first 100 of {len(df)} rows)")
            else:
                log.info(
                    "\n" + tabulate.tabulate(df, headers="keys", tablefmt="simple")
                )
        else:
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
        log.info(f"Successfully wrote {len(df)} sequences to {output}")
