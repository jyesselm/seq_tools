"""
commandline interface for seq_tools
"""

import os
import click
import tabulate
import pandas as pd

from seq_tools import sequence, dataframe
from seq_tools.logger import setup_applevel_logger, get_logger
from seq_tools.validation import validate_dataframe as validate_df, ensure_name_column
from seq_tools.utils import get_resources_path

pd.set_option("display.max_colwidth", None)

log = get_logger("cli")


def validate_dataframe(df) -> None:
    """
    validates a dataframe to have a column named `sequence` and `name`
    :param df: dataframe with sequences
    :return: None
    """
    validate_df(df, require_name=False)
    ensure_name_column(df)


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


def handle_output(df, output, show_all: bool = False) -> None:
    """
    handles the output of the dataframe
    :param df: dataframe with sequences
    :param output: output file
    :param show_all: if True, show all sequences; if False, show only first
    :return: None
    """
    log = get_logger("handle_output")
    if len(df) == 1:
        log.info(f"output->\n{df.iloc[0]}")
        log.info(f"Writing output to: {output}")
        df.to_csv(output, index=False)
    else:
        log.info(f"Writing output CSV to: {output}")
        if show_all:
            if len(df) > 100:
                log.info(
                    "\n" + tabulate.tabulate(df[0:100], headers="keys", tablefmt="simple")
                )
                log.info(f"... (showing first 100 of {len(df)} rows)")
            else:
                log.info("\n" + tabulate.tabulate(df, headers="keys", tablefmt="simple"))
        else:
            log.info(f"\nFirst sequence (showing 1 of {len(df)}):")
            log.info("\n" + tabulate.tabulate(
                df.iloc[0:1],
                headers="keys",
                tablefmt="simple",
                showindex=False
            ))
        df.to_csv(output, index=False)
        log.info(f"Successfully wrote {len(df)} sequences to {output}")


@click.group()
def cli():
    """
    a set scripts to manipulate sequences in csv files
    """


@cli.command(help="add a sequence to 5' and/or 3'")
@click.argument("data")
@click.option("-p5", "--p5-seq", default="")
@click.option("-p3", "--p3-seq", default="")
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def add(data, p5_seq, p3_seq, output):
    """
    adds a sequence to a dataframe
    :param data: can be a sequence or a file
    :param p5_seq: sequence to add to 5'
    :param p3_seq: sequence to add to 3'
    :param output: output file
    """
    setup_applevel_logger()
    log = get_logger("add")
    df = get_input_dataframe(data)
    
    if p5_seq:
        log.info(f"Adding 5' sequence: {p5_seq} (length: {len(p5_seq)})")
    if p3_seq:
        log.info(f"Adding 3' sequence: {p3_seq} (length: {len(p3_seq)})")
    
    df = dataframe.add(df, p5_seq, p3_seq)
    handle_output(df, output)


@cli.command(help="calculate the edit distance of a library")
@click.option(
    "-p", "--parallel", is_flag=True, help="calculate the edit distance in parallel"
)
@click.option("-w", "--workers", type=int, help="number of workers to use", default=1)
@click.option(
    "--use-threads",
    is_flag=True,
    help="use threads instead of processes",
    default=False,
)
@click.argument("data", type=click.Path(exists=True))
def edit_distance(data, parallel, workers, use_threads):
    """
    calculates the edit distance of a library
    :param data: can be a sequence or a file
    """
    setup_applevel_logger()
    df = pd.read_csv(data)
    if parallel:
        log.info(f"calculating edit distance in parallel with {workers} workers")
        if use_threads:
            log.info("using threads instead of processes")
        score = dataframe.calc_edit_distance_parallel(df, workers, use_threads)
    else:
        score = dataframe.calc_edit_distance(df)
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
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
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
    if len(df) > 0 and "extinction_coeff" in df.columns:
        first_ec = df.iloc[0]["extinction_coeff"]
        log.info(f"Calculated extinction coefficients for {ntype} ({'double-stranded' if double_stranded else 'single-stranded'})")
        if len(df) == 1:
            log.info(f"Extinction coefficient: {first_ec:.0f} M⁻¹cm⁻¹")
        else:
            log.info(f"First sequence extinction coefficient: {first_ec:.0f} M⁻¹cm⁻¹")
            log.info(f"Average extinction coefficient: {df['extinction_coeff'].mean():.0f} M⁻¹cm⁻¹")
    handle_output(df, output)


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
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def mw(data, ntype, double_stranded, output):
    """
    calculates the molecular weight for each sequence
    :param data:
    :param double_stranded:
    :param output:
    :return:
    """
    setup_applevel_logger()
    log = get_logger("molecular_weight")
    df = get_input_dataframe(data)
    ntype = get_ntype(df, ntype)
    df = dataframe.get_molecular_weight(df, ntype, double_stranded)
    if len(df) > 0 and "mw" in df.columns:
        first_mw = df.iloc[0]["mw"]
        log.info(f"Calculated molecular weights for {ntype} ({'double-stranded' if double_stranded else 'single-stranded'})")
        if len(df) == 1:
            log.info(f"Molecular weight: {first_mw:.2f} Da")
        else:
            log.info(f"First sequence molecular weight: {first_mw:.2f} Da")
            log.info(f"Average molecular weight: {df['mw'].mean():.2f} Da")
    handle_output(df, output)


@cli.command(help="calculate reverse complement for each sequence")
@click.argument("data")
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def rc(data, ntype, output):
    """
    calculates the reverse complement for each sequence
    :param data: can be a sequence or a file
    :param output: output file
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


@cli.command(help="fold rna sequences")
@click.argument("data")
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def fold(data, output):
    """
    fold rna sequences
    :param data: can be a sequence or a file
    """
    setup_applevel_logger()
    log = get_logger("fold")
    df = get_input_dataframe(data)
    df = dataframe.fold(df)
    if len(df) > 0 and "structure" in df.columns:
        first_row = df.iloc[0]
        log.info(f"Folded sequences using ViennaRNA")
        if len(df) == 1:
            log.info(f"Structure: {first_row.get('structure', 'N/A')}")
            log.info(f"MFE: {first_row.get('mfe', 'N/A')} kcal/mol")
            log.info(f"Ensemble defect: {first_row.get('ens_defect', 'N/A')}")
    handle_output(df, output)


@cli.command(help="checks to see if p5 is present in all sequences")
@click.argument("data")
@click.option("-p5", "--p5-seq", help="p5 sequence (if not provided, finds longest common prefix)", default=None)
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
def has_p5(data, p5_seq, ntype):
    """
    checks if a sequence has a p5 sequence. If no p5 sequence is provided,
    automatically finds the longest common prefix from all sequences.
    :param data: can be a sequence or a file
    :param p5_seq: p5 sequence (optional)
    :param ntype: type of nucleic acid
    """
    setup_applevel_logger()
    log = get_logger("has_p5")
    df = get_input_dataframe(data)
    if ntype is not None:
        get_ntype(df, ntype)
    
    # Load known p5 sequences
    df_p5 = pd.read_csv(get_resources_path() / "p5_sequences.csv")
    
    if p5_seq is None:
        # Find longest common prefix
        common_prefix = dataframe.find_longest_common_prefix(df)
        if common_prefix:
            log.info(f"Found longest common 5' sequence: {common_prefix} (length: {len(common_prefix)})")
            log.info("This sequence is present at the 5' end of all sequences")
            
            # Check if it matches any known p5 sequences
            matching_known = []
            for _, row in df_p5.iterrows():
                known_seq = row["sequence"]
                known_name = row.get("name", "unknown")
                # Check if the common prefix starts with or matches a known sequence
                if common_prefix.startswith(known_seq):
                    matching_known.append(row)
                elif known_seq.startswith(common_prefix):
                    matching_known.append(row)
            
            if matching_known:
                log.info("\nMatching known p5 sequences:")
                for row in matching_known:
                    name = row.get("name", "unknown")
                    seq = row.get("sequence", "")
                    code = row.get("code", "")
                    structure = row.get("structure", "")
                    log.info(f"  - {name}: {seq}")
                    if code:
                        log.info(f"    Code: {code}")
                    if structure:
                        log.info(f"    Structure: {structure}")
            else:
                log.info("No matching known p5 sequences found in resource file")
        else:
            log.info("No common 5' sequence found in all sequences")
    else:
        has_p5_seq = dataframe.has_5p_sequence(df, p5_seq)
        if has_p5_seq:
            log.info(f"p5 sequence '{p5_seq}' is present in all sequences")
            
            # Check if it matches any known p5 sequences
            matching_known = []
            for _, row in df_p5.iterrows():
                known_seq = row["sequence"]
                if p5_seq == known_seq or p5_seq.startswith(known_seq) or known_seq.startswith(p5_seq):
                    matching_known.append(row)
            
            if matching_known:
                log.info("\nMatching known p5 sequences:")
                for row in matching_known:
                    name = row.get("name", "unknown")
                    seq = row.get("sequence", "")
                    code = row.get("code", "")
                    structure = row.get("structure", "")
                    log.info(f"  - {name}: {seq}")
                    if code:
                        log.info(f"    Code: {code}")
                    if structure:
                        log.info(f"    Structure: {structure}")
        else:
            log.info(f"p5 sequence '{p5_seq}' is not present in all sequences")


@cli.command(help="checks to see if p3 is present in all sequences")
@click.argument("data")
@click.option("-p3", "--p3-seq", help="p3 sequence (if not provided, finds longest common suffix)", default=None)
@click.option(
    "-nt",
    "--ntype",
    default=None,
    type=click.Choice([None, "RNA", "DNA"]),
    help="type of nucleic acid",
)
def has_p3(data, p3_seq, ntype):
    """
    checks if a sequence has a p3 sequence. If no p3 sequence is provided,
    automatically finds the longest common suffix from all sequences.
    :param data: can be a sequence or a file
    :param p3_seq: p3 sequence (optional)
    :param ntype: type of nucleic acid
    """
    setup_applevel_logger()
    log = get_logger("has_p3")
    df = get_input_dataframe(data)
    if ntype is not None:
        get_ntype(df, ntype)
    
    # Load known p3 sequences
    df_p3 = pd.read_csv(get_resources_path() / "p3_sequences.csv")
    
    if p3_seq is None:
        # Find longest common suffix
        common_suffix = dataframe.find_longest_common_suffix(df)
        if common_suffix:
            log.info(f"Found longest common 3' sequence: {common_suffix} (length: {len(common_suffix)})")
            log.info("This sequence is present at the 3' end of all sequences")
            
            # Check if it matches any known p3 sequences
            matching_known = []
            for _, row in df_p3.iterrows():
                known_seq = row["sequence"]
                # Check if the common suffix ends with or matches a known sequence
                if common_suffix.endswith(known_seq):
                    matching_known.append(row)
                elif known_seq.endswith(common_suffix):
                    matching_known.append(row)
            
            if matching_known:
                log.info("\nMatching known p3 sequences:")
                for row in matching_known:
                    name = row.get("name", "unknown")
                    seq = row.get("sequence", "")
                    code = row.get("code", "")
                    structure = row.get("structure", "")
                    log.info(f"  - {name}: {seq}")
                    if code:
                        log.info(f"    Code: {code}")
                    if structure:
                        log.info(f"    Structure: {structure}")
            else:
                log.info("No matching known p3 sequences found in resource file")
        else:
            log.info("No common 3' sequence found in all sequences")
    else:
        has_p3_seq = dataframe.has_3p_sequence(df, p3_seq)
        if has_p3_seq:
            log.info(f"p3 sequence '{p3_seq}' is present in all sequences")
            
            # Check if it matches any known p3 sequences
            matching_known = []
            for _, row in df_p3.iterrows():
                known_seq = row["sequence"]
                if p3_seq == known_seq or p3_seq.endswith(known_seq) or known_seq.endswith(p3_seq):
                    matching_known.append(row)
            
            if matching_known:
                log.info("\nMatching known p3 sequences:")
                for row in matching_known:
                    name = row.get("name", "unknown")
                    seq = row.get("sequence", "")
                    code = row.get("code", "")
                    structure = row.get("structure", "")
                    log.info(f"  - {name}: {seq}")
                    if code:
                        log.info(f"    Code: {code}")
                    if structure:
                        log.info(f"    Structure: {structure}")
        else:
            log.info(f"p3 sequence '{p3_seq}' is not present in all sequences")


@cli.command(help="convert rna sequence(s) to dna")
@click.argument("data")
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def to_dna(data, output):
    """
    Convert RNA sequence to DNA
    """
    setup_applevel_logger()
    log = get_logger("to_dna")
    df = get_input_dataframe(data)
    original_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
    df = df[["name", "sequence"]]
    df = dataframe.to_dna(df)
    if original_seq:
        converted_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
        log.info(f"Converted RNA to DNA (U -> T)")
        if len(df) == 1:
            log.info(f"Example: {original_seq[:20]}... -> {converted_seq[:20]}...")
    handle_output(df, output)


@cli.command(help="convert rna sequence(s) to dna template, includes T7 promoter")
@click.argument("data")
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def to_dna_template(data, output):
    """
    Convert RNA sequence to DNA template with T7 promoter
    """
    setup_applevel_logger()
    log = get_logger("to_dna_template")
    df = get_input_dataframe(data)
    original_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
    df = df[["name", "sequence"]]
    df = dataframe.to_dna_template(df)
    if original_seq:
        converted_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
        log.info(f"Converted RNA to DNA template with T7 promoter")
        if len(df) == 1:
            log.info(f"Added T7 promoter (20 bases) and converted U -> T")
            log.info(f"Example: {original_seq[:20]}... -> {converted_seq[:40]}...")
    handle_output(df, output)


@cli.command(help="generate fasta file from csv")
@click.argument("data")
@click.option("-o", "--output", help="output file (default: output.fasta)", default="output.fasta")
def to_fasta(data, output):
    """
    generate fasta file from csv
    :param data: can be a sequence or a file
    :param output: output file
    """
    setup_applevel_logger()
    log = get_logger("to_fasta")
    df = get_input_dataframe(data)
    log.info(f"Writing FASTA output to: {output}")
    dataframe.to_fasta(df, output)
    log.info(f"Successfully wrote {len(df)} sequences to {output}")


@cli.command(help="generate oligo pool file from csv")
@click.argument("data")
@click.option("-n", "--name", help="name of the opool file", default="opool")
@click.option("-o", "--output", help="output file (default: output.xlsx)", default="output.xlsx")
def to_opool(data, name, output):
    """
    generate opool file from csv
    :param data: can be a sequence or a file
    :param output: output file
    """
    setup_applevel_logger()
    log = get_logger("to_opool")
    df = get_input_dataframe(data)
    log.info(f"Writing opool output to: {output}")
    dataframe.to_opool(df, name, output)
    log.info(f"Successfully wrote {len(df)} sequences to {output}")


@cli.command(help="convert dna sequence(s) to rna")
@click.argument("data")
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def to_rna(data, output):
    """
    Convert DNA sequence to RNA
    """
    setup_applevel_logger()
    log = get_logger("to_rna")
    df = get_input_dataframe(data)
    original_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
    df = df[["name", "sequence"]]
    # apply sequence.to_rna to `sequence` column
    df["sequence"] = df["sequence"].apply(sequence.to_rna)
    if original_seq:
        converted_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
        log.info(f"Converted DNA to RNA (T -> U)")
        if len(df) == 1:
            log.info(f"Example: {original_seq[:20]}... -> {converted_seq[:20]}...")
    handle_output(df, output)


@cli.command(help="trim 5'/3' ends of sequences")
@click.argument("data")
@click.option("-p5", "--p5-cut", default=0)
@click.option("-p3", "--p3-cut", default=0)
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def trim(data, p5_cut, p3_cut, output):
    """
    trim 5'/3' ends of sequences
    :param data: can be a sequence or a file
    :param p5_cut: trim off 5' end
    :param p3_cut: trim off 3' end
    :param output: output file
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


@cli.command(help="convert dna sequence(s) to rna")
@click.argument("data")
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def transcribe(data, output):
    """
    Transcribe DNA sequence to RNA (removes T7 promoter if present)
    """
    setup_applevel_logger()
    log = get_logger("transcribe")
    df = get_input_dataframe(data)
    original_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
    df = df[["name", "sequence"]]
    df = dataframe.transcribe(df)
    if original_seq:
        transcribed_seq = df.iloc[0]["sequence"] if len(df) > 0 else ""
        log.info(f"Transcribed DNA to RNA (removed T7 promoter if present, T -> U)")
        if len(df) == 1:
            log.info(f"Example: {original_seq[:30]}... -> {transcribed_seq[:20]}...")
    handle_output(df, output)


@cli.command(
    help="generate mutated sequences from a template with optional constant 5' and 3' ends"
)
@click.argument("template")
@click.option(
    "-n",
    "--num-sequences",
    type=int,
    default=10,
    help="number of mutated sequences to generate",
)
@click.option(
    "-m",
    "--num-mutations",
    type=int,
    required=True,
    help="number of mutations per sequence",
)
@click.option(
    "-p5",
    "--p5-seq",
    default=None,
    help="constant 5' sequence (mutations only in middle region)",
)
@click.option(
    "-p3",
    "--p3-seq",
    default=None,
    help="constant 3' sequence (mutations only in middle region)",
)
@click.option(
    "-nt",
    "--ntype",
    default="DNA",
    type=click.Choice(["DNA", "RNA"]),
    help="type of nucleic acid",
)
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def mutate(template, num_sequences, num_mutations, p5_seq, p3_seq, ntype, output):
    """
    generates mutated sequences from a template sequence with optional constant 5' and 3' ends

    :param template: template sequence to mutate (can be a sequence string or file with one sequence)
    :param num_sequences: number of mutated sequences to generate
    :param num_mutations: number of mutations to introduce per sequence
    :param p5_seq: optional constant 5' sequence (if provided, mutations only in middle region)
    :param p3_seq: optional constant 3' sequence (if provided, mutations only in middle region)
    :param ntype: type of nucleic acid (DNA or RNA)
    :param output: output file
    """
    setup_applevel_logger()
    log = get_logger("mutate")

    # Get template sequence - can be a file or a sequence string
    if os.path.isfile(template):
        log.info(f"reading template from file {template}")
        df_temp = pd.read_csv(template)
        validate_dataframe(df_temp)
        if len(df_temp) != 1:
            raise ValueError(
                f"Template file must contain exactly one sequence, found {len(df_temp)}"
            )
        template_seq = df_temp.iloc[0]["sequence"]
    else:
        log.info(f"using template sequence: {template}")
        template_seq = template

    log.info(
        f"generating {num_sequences} sequences with {num_mutations} mutations each"
    )
    if p5_seq:
        log.info(f"using constant 5' sequence: {p5_seq}")
    if p3_seq:
        log.info(f"using constant 3' sequence: {p3_seq}")

    df = dataframe.generate_mutated_sequences(
        template=template_seq,
        num_mutations=num_mutations,
        num_sequences=num_sequences,
        p5_seq=p5_seq,
        p3_seq=p3_seq,
        ntype=ntype,
    )
    handle_output(df, output)


@cli.command(help="generate random sequences with optional constant 5' and 3' ends")
@click.option(
    "-l",
    "--length",
    type=int,
    required=True,
    help="total length of sequences (including constant 5' and 3' if provided)",
)
@click.option(
    "-n",
    "--num-sequences",
    type=int,
    default=10,
    help="number of random sequences to generate",
)
@click.option(
    "-p5",
    "--p5-seq",
    default=None,
    help="constant 5' sequence",
)
@click.option(
    "-p3",
    "--p3-seq",
    default=None,
    help="constant 3' sequence",
)
@click.option(
    "-nt",
    "--ntype",
    default="DNA",
    type=click.Choice(["DNA", "RNA"]),
    help="type of nucleic acid",
)
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def random(length, num_sequences, p5_seq, p3_seq, ntype, output):
    """
    generates random sequences with optional constant 5' and 3' ends

    :param length: total length of sequences (including constant 5' and 3' if provided)
    :param num_sequences: number of random sequences to generate
    :param p5_seq: optional constant 5' sequence
    :param p3_seq: optional constant 3' sequence
    :param ntype: type of nucleic acid (DNA or RNA)
    :param output: output file
    """
    setup_applevel_logger()
    log = get_logger("random")

    log.info(f"generating {num_sequences} random sequences of length {length}")
    if p5_seq:
        log.info(f"using constant 5' sequence: {p5_seq}")
    if p3_seq:
        log.info(f"using constant 3' sequence: {p3_seq}")

    df = dataframe.generate_random_sequences(
        length=length,
        num_sequences=num_sequences,
        p5_seq=p5_seq,
        p3_seq=p3_seq,
        ntype=ntype,
    )
    handle_output(df, output)


@cli.command(name="list-common-seqs", help="list available 5' and 3' common sequences")
def list_common_seqs():
    """
    Lists all available 5' and 3' common sequences from resource files in nice tables.
    """
    setup_applevel_logger()
    log = get_logger("list_common_seqs")
    
    # Load p5 sequences
    p5_path = get_resources_path() / "p5_sequences.csv"
    if p5_path.exists():
        df_p5 = pd.read_csv(p5_path)
        log.info("\nAvailable 5' sequences:")
        log.info("\n" + tabulate.tabulate(
            df_p5,
            headers="keys",
            tablefmt="simple",
            showindex=False
        ))
    else:
        log.warning(f"p5_sequences.csv not found at {p5_path}")
    
    # Load p3 sequences
    p3_path = get_resources_path() / "p3_sequences.csv"
    if p3_path.exists():
        df_p3 = pd.read_csv(p3_path)
        log.info("\n\nAvailable 3' sequences:")
        log.info("\n" + tabulate.tabulate(
            df_p3,
            headers="keys",
            tablefmt="simple",
            showindex=False
        ))
    else:
        log.warning(f"p3_sequences.csv not found at {p3_path}")


@cli.command(name="remove-common-seqs", help="identify and remove common 5' and 3' sequences, checking both sequence and structure")
@click.argument("data")
@click.option("-o", "--output", help="output file (default: output.csv)", default="output.csv")
def remove_common_seqs(data, output):
    """
    Identifies and removes common 5' and 3' sequences from sequences in a CSV file.
    
    This command checks all sequences in the input CSV against known 5' and 3'
    sequences from resource files. It attempts to match both sequence and structure
    patterns, and provides warnings if only sequence or only structure matches.
    
    :param data: Input CSV file with sequences and optionally structures
    :param output: Output CSV file
    """
    setup_applevel_logger()
    log = get_logger("remove_common_seqs")
    
    df = get_input_dataframe(data)
    
    try:
        df, result = dataframe.remove_common_seqs(df)
        
        # Display which sequences were removed
        removed = result.get("removed", {})
        if removed.get("p5_sequence"):
            log.info(f"Removed 5' sequence: {removed['p5_sequence']} (length: {removed['p5_length']})")
        if removed.get("p3_sequence"):
            log.info(f"Removed 3' sequence: {removed['p3_sequence']} (length: {removed['p3_length']})")
        
        # Display warnings if any
        if result.get("warnings"):
            for warning in result["warnings"]:
                log.warning(warning)
        
        log.info("Successfully identified and removed common 5' and/or 3' sequences")
        
        # Show only first sequence in output
        if len(df) > 0:
            log.info(f"\nFirst sequence (showing 1 of {len(df)}):")
            log.info("\n" + tabulate.tabulate(
                df.iloc[0:1],
                headers="keys",
                tablefmt="simple",
                showindex=False
            ))
        df.to_csv(output, index=False)
        log.info(f"\nSuccessfully wrote {len(df)} sequences to {output}")
    except ValueError as e:
        log.error(f"Error: {e}")
        raise click.ClickException(str(e))


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    cli()
