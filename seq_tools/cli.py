"""
commandline interface for seq_tools
"""
import os
import click
import tabulate
import pandas as pd

from seq_tools import sequence, dataframe
from seq_tools.logger import setup_applevel_logger, get_logger

pd.set_option("display.max_colwidth", None)


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
    log = get_logger("get_input_dataframe")
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
    log = get_logger("handle_output")
    if len(df) == 1:
        log.info(f"output->\n{df.iloc[0]}")
    else:
        log.info(f"output csv: {output}")
        if len(df) > 100:
            log.info(
                "\n"
                + tabulate.tabulate(
                    df[0:100], headers="keys", tablefmt="simple"
                )
            )
        else:
            log.info(
                "\n" + tabulate.tabulate(df, headers="keys", tablefmt="simple")
            )
        df.to_csv(output, index=False)


@click.group()
def cli():
    """
    main function for script
    """

@cli.command(help="add a sequence to 5' and/or 3'")
@click.argument("data")
@click.option("-p5", "--p5-seq", default="")
@click.option("-p3", "--p3-seq", default="")
@click.option("-o", "--output", help="output file", default="output.csv")
def add(data, p5_seq, p3_seq, output):
    """
    adds a sequence to a dataframe
    :param data: can be a sequence or a file
    :param p5_seq: sequence to add to 5'
    :param p3_seq: sequence to add to 3'
    :param output: output file
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    df = dataframe.add(df, p5_seq, p3_seq)
    handle_output(df, output)

@cli.command(help="fold rna sequences")
@click.argument("data")
@click.option("-o", "--output", help="output file", default="output.csv")
def fold(data, output):
    """
    fold rna sequences
    :param data: can be a sequence or a file
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    df = dataframe.fold(df)
    handle_output(df, output)


@cli.command(help="convert rna sequence(s) to dna")
@click.argument("data")
@click.option("-o", "--output", help="output file", default="output.csv")
def to_dna(data, output):
    """
    Convert RNA sequence to DNA
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    df = df[["name", "sequence"]]
    df = dataframe.to_dna(df)
    handle_output(df, output)


@cli.command(
    help="convert rna sequence(s) to dna template, includes T7 promoter"
)
@click.argument("data")
@click.option("-o", "--output", help="output file", default="output.csv")
def to_dna_template(data, output):
    """
    Convert RNA sequence to DNA
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    df = df[["name", "sequence"]]
    df = dataframe.to_dna(df)
    handle_output(df, output)


@cli.command(help="convert rna sequence(s) to dna")
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
    df["sequence"] = df["sequence"].apply(sequence.to_rna)
    handle_output(df, output)


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    cli()
