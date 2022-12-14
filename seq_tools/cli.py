"""
commandline interface for seq_tools
"""
import os
import click
import pandas as pd

from seq_tools import sequence
from seq_tools.logger import setup_applevel_logger, get_logger

pd.set_option('display.max_colwidth', None)


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
        log.info(f"csv file contains {len(df)}")
    else:
        log.info(f"reading sequence {data}")
        data_df = [["seq", data]]
        df = pd.DataFrame(data_df, columns=["name", "sequence"])
    validate_dataframe(df)
    return df


def handle_output(df, output) -> None:
    """
    handles the output of the dataframe
    :param df: dataframe with sequences
    :param output: output file
    :return: None
    """
    log = get_logger('handle_output')
    if len(df) == 1:
        log.info(f"output->\n{df.iloc[0]}")
    else:
        log.info(f"output csv: {output}")
        df.to_csv(output, index=False)


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
    df = get_input_dataframe(data)
    df = df[["name", "sequence"]]
    # apply sequence.to_dna to `sequence` column
    df['sequence'] = df['sequence'].apply(sequence.to_dna)
    handle_output(df, output)


@cli.command(help='convert rna sequence(s) to dna')
@click.argument("data")
@click.option("-o", "--output", help="output file", default="output.csv")
def to_rna(data, output):
    """
    Convert DNA sequence to RNA
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    df = df[["name", "sequence"]]
    # apply sequence.to_dna to `sequence` column
    df['sequence'] = df['sequence'].apply(sequence.to_rna)
    handle_output(df, output)


# pylint: disable=no-value-for-parameter
if __name__ == '__main__':
    cli()
