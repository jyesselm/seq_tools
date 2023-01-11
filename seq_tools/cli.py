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
        log.info(f"csv file contains {len(df)} sequences")
    else:
        log.info(f"reading sequence {data}")
        data_df = [["seq", data]]
        df = pd.DataFrame(data_df, columns=["name", "sequence"])
    validate_dataframe(df)
    return df


def get_ntype(df, ntype) -> str:
    """
    handles the ntype parameter
    :param df: dataframe with sequences
    :param ntype: nucleotide type
    :return: str
    """
    log = get_logger("handle_ntype")
    df_ntype = dataframe.determine_ntype(df)
    log.info(f"determining nucleic acid type: {df_ntype}")
    # enforce ntype
    if ntype == "DNA":
        log.info("forcing sequences to be DNA")
        dataframe.to_dna(df)
        return ntype
    if ntype == "RNA":
        log.info("forcing sequences to be RNA")
        dataframe.to_rna(df)
        return ntype
    return df_ntype


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
    a set scripts to manipulate sequences in csv files
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


@cli.command(help="calculate the edit distance of a library")
@click.argument("data", type=click.Path(exists=True))
def edit_distance(data):
    """
    calculates the edit distance of a library
    :param data: can be a sequence or a file
    """
    setup_applevel_logger()
    df = pd.read_csv(data)
    score = dataframe.calc_edit_distance(df)
    log = get_logger("edit_distance")
    log.info(f"edit distance: {score}")


@cli.command(help="calculate the extinction coefficient for each sequence")
@click.argument("data")
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
@click.option("-ds", "--double-stranded", is_flag=True)
@click.option("-o", "--output", help="output file", default="output.csv")
def ec(data, ntype, double_stranded, output):
    """
    calculates the extinction coefficient for each sequence
    :param data: can be a sequence or a file
    :param ntype: type of nucleic acid
    :param double_stranded: if the sequence is double stranded
    :param output: output file
    """
    setup_applevel_logger()
    log = get_logger("extinction_coeff")
    df = get_input_dataframe(data)
    ntype = get_ntype(df, ntype)
    df = dataframe.get_extinction_coeff(df, ntype, double_stranded)
    handle_output(df, output)
    if len(df) != 1:
        log.info(
            "avg extinction coefficient: " + str(df["extinction_coeff"].mean())
        )


@cli.command(help="calculate the molecular weight for each sequence")
@click.argument("data")
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
@click.option("-ds", "--double-stranded", is_flag=True)
@click.option("-o", "--output", help="output file", default="output.csv")
def mw(data, ntype, double_stranded, output):
    """
    calculates the molecular weight for each sequence
    :param data:
    :param double_stranded:
    :param output:
    :return:
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    ntype = get_ntype(df, ntype)
    df = dataframe.get_molecular_weight(df, ntype, double_stranded)
    handle_output(df, output)
    log = get_logger("molecular_weight")
    if len(df) != 1:
        log.info("avg molecular weight: " + str(df["molecular_weight"].mean()))


@cli.command(help="calculate reverse complement for each sequence")
@click.argument("data")
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
@click.option("-o", "--output", help="output file", default="output.csv")
def rc(data, ntype, output):
    """
    calculates the reverse complement for each sequence
    :param data: can be a sequence or a file
    :param output: output file
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    ntype = get_ntype(df, ntype)
    df = dataframe.get_reverse_complement(df, ntype)
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


@cli.command(help="checks to see if p5 is present in all sequences")
@click.argument("data")
@click.option("-p5", "--p5-seq", help="p5 sequence", required=True)
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
def has_p5(data, p5_seq, ntype):
    """
    checks if a sequence has a p5 sequence
    :param data: can be a sequence or a file
    :param p5_seq: p5 sequence
    :param ntype: type of nucleic acid
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    get_ntype(df, ntype)
    has_p5_seq = dataframe.has_5p_sequence(df, p5_seq)
    log = get_logger("has_p5")
    if has_p5_seq:
        log.info("p5 sequence is present in all sequences")
    else:
        log.info("p5 sequence is not present in all sequences")


@cli.command(help="checks to see if p3 is present in all sequences")
@click.argument("data")
@click.option("-p3", "--p3-seq", help="p3 sequence", required=True)
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
def has_p3(data, p3_seq, ntype):
    """
    checks if a sequence has a p3 sequence
    :param data: can be a sequence or a file
    :param p3_seq: p3 sequence
    :param ntype: type of nucleic acid
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    get_ntype(df, ntype)
    has_p3_seq = dataframe.has_5p_sequence(df, p3_seq)
    log = get_logger("has_p3")
    if has_p3_seq:
        log.info("p3 sequence is present in all sequences")
    else:
        log.info("p3 sequence is not present in all sequences")


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
    df = dataframe.to_dna_template(df)
    handle_output(df, output)


@cli.command(help="generate fasta file from csv")
@click.argument("data")
@click.option("-o", "--output", help="output file", default="test.fasta")
def to_fasta(data, output):
    """
    generate fasta file from csv
    :param data: can be a sequence or a file
    :param output: output file
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    dataframe.to_fasta(df, output)


@cli.command(help="generate oligo pool file from csv")
@click.argument("data")
@click.option("-n", "--name", help="name of the opool file", default="opool")
@click.option("-o", "--output", help="output file", default="opool.xlsx")
def to_opool(data, name, output):
    """
    generate opool file from csv
    :param data: can be a sequence or a file
    :param output: output file
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    dataframe.to_opool(df, name, output)


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


@cli.command(help="trim 5'/3' ends of sequences")
@click.argument("data")
@click.option("-p5", "--p5-cut", default=0)
@click.option("-p3", "--p3-cut", default=0)
@click.option("-o", "--output", help="output file", default="output.csv")
def trim(data, p5_cut, p3_cut, output):
    """
    trim 5'/3' ends of sequences
    :param data: can be a sequence or a file
    :param p5_cut: trim off 5' end
    :param p3_cut: trim off 3' end
    :param output: output file
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    df = dataframe.trim(df, p5_cut, p3_cut)
    handle_output(df, output)


@cli.command(help="convert dna sequence(s) to rna")
@click.argument("data")
@click.option("-o", "--output", help="output file", default="output.csv")
def transcribe(data, output):
    """
    Convert DNA sequence to RN
    """
    setup_applevel_logger()
    df = get_input_dataframe(data)
    df = df[["name", "sequence"]]
    df = dataframe.transcribe(df)
    handle_output(df, output)


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    cli()
