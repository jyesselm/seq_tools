"""
commandline interface for seq_tools
"""
import os
import click
import pandas as pd

from seq_tools import sequence
from seq_tools.logger import setup_applevel_logger, get_logger


def validate_dataframe(df) -> None:
    """
    validates a dataframe to have a column named `sequence` and `name`
    :param df: dataframe with sequences
    :return: None
    """
    if "sequence" not in df.columns:
        raise ValueError("sequence column not found")
    if "name" not in df.columns:
        df["name"] = [f"seq_{i}" for i in range(len(df))]


def get_input_dataframe(data) -> pd.DataFrame:
    """
    returns a dataframe from a sequence or a file
    :param data: can be a seqeunce or a file
    :return: pd.DataFrame
    """
    log = get_logger('get_input_dataframe')
    if os.path.isfile(data):
        log.info(f"reading file {data}")
        df = pd.read_csv(data)
    else:
        log.info(f"reading sequence {data}")
        data_df = [["seq", data]]
        df = pd.DataFrame(data_df, columns=["name", "sequence"])
    validate_dataframe(df)
    return df


@click.group()
def cli():
    """
    main function for script
    """


@cli.command(help='convert rna sequence(s) to dna')
@click.argument("data")
@click.option("-o", "--output", help="output file", default="output.csv")
def to_dna(data, output):
    """
    Convert RNA sequence to DNA
    """
    setup_applevel_logger()
    log = get_logger('to_dna')
    df = get_input_dataframe(data)
    # apply sequence.to_dna to `sequence` column
    df['sequence'] = df['sequence'].apply(sequence.to_dna)
    if len(df) == 1:
        log.info(f"converted sequence: {df.iloc[0]['sequence']}")
    else:
        log.info(f"converted sequences: {output}")
        df.to_csv(output, index=False)


# pylint: disable=no-value-for-parameter
if __name__ == '__main__':
    cli()
