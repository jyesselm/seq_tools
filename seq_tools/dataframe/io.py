"""I/O operations for sequence dataframes."""

import pandas as pd


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
    df.to_excel(filename, index=False)
