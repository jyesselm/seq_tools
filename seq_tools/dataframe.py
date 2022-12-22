"""
module for working with dataframes that contain nucleotide sequences
"""
import pandas as pd
import numpy as np
import editdistance
import vienna

from seq_tools import sequence, extinction_coeff


def add(df: pd.DataFrame, p5_seq: str, p3_seq: str) -> pd.DataFrame:
    """
    adds a 5' and 3' sequence to the sequences in the dataframe
    :param df: dataframe
    :param p5_seq: 5' sequence
    :param p3_seq: 3' sequence
    :return: None
    """
    df["sequence"] = df["sequence"].apply(lambda x: p5_seq + x + p3_seq)
    if "structure" in df.columns:
        df = fold(df)
    return df


def calc_edit_distance(df: pd.DataFrame) -> float:
    """
    calculates the edit distance between each sequence in the dataframe
    :param df: dataframe
    :return: the edit distance
    """
    if len(df) == 1:
        return 0
    scores = [100 for _ in range(len(df))]
    sequences = list(df["sequence"])
    for i, seq1 in enumerate(sequences):
        for j, seq2 in enumerate(sequences):
            if i >= j:
                continue
            diff = editdistance.eval(seq1, seq2)
            if scores[i] > diff:
                scores[i] = diff
            if scores[j] > diff:
                scores[j] = diff
    avg = np.mean(scores)
    return avg


def determine_ntype(df: pd.DataFrame) -> str:
    """
    determines the nucleotide type of the sequences in the dataframe
    :param df: dataframe
    :return: nucleotide type, RNA or DNA
    """
    results = []
    for _, row in df.iterrows():
        ntype = "UNCERTAIN"
        if row["sequence"].count("T") > 0:
            ntype = "DNA"
        elif row["sequence"].count("U") > 0:
            ntype = "RNA"
        results.append(ntype)
    if df["sequence"].str.len().mean() > 10:
        if results.count("DNA") > 0 and results.count("RNA") > 0:
            raise ValueError("Cannot determine nucleotide type")
    if results.count("RNA") > 0:
        return "RNA"
    return "DNA"


def fold(df: pd.DataFrame) -> pd.DataFrame:
    """
    folds each sequence in the dataframe
    :param df: dataframe
    """

    def _fold(seq):
        v_res = vienna.fold(seq)
        return pd.Series(
            [v_res.dot_bracket, v_res.mfe, v_res.ens_defect],
            index=["structure", "mfe", "ens_defect"],
        )

    df[["structure", "mfe", "ens_defect"]] = df["sequence"].apply(_fold)
    return df


def get_extinction_coeff(
    df: pd.DataFrame, ntype: str, double_stranded: bool
) -> pd.DataFrame:
    """
    calculates the extinction coefficient for each sequence in the dataframe
    :param df: dataframe
    :param ntype: nucleotide type, RNA or DNA
    :param double_stranded: is double stranded?
    :return: None
    """

    def compute_w_struct(row) -> float:
        """
        computes the extinction coefficient for a sequence with a structure
        :param row: dataframe row
        :return: extinction coefficient
        """
        return extinction_coeff.get_extinction_coeff(
            row["sequence"], ntype, double_stranded, row["structure"]
        )

    if ntype == "RNA" and "structure" in df.columns:
        df["extinction_coeff"] = df.apply(compute_w_struct, axis=1)
    else:
        df["extinction_coeff"] = df["sequence"].apply(
            lambda x: extinction_coeff.get_extinction_coeff(
                x, ntype, double_stranded
            )
        )
    return df


def get_length(df: pd.DataFrame) -> pd.DataFrame:
    """
    calculates the length of each sequence in the dataframe
    :param df: dataframe
    :return: None
    """
    df["length"] = df["sequence"].apply(len)
    return df


def get_molecular_weight(
    df: pd.DataFrame, ntype: str, double_stranded: bool
) -> pd.DataFrame:
    """
    :param df: pandas data frame
    :param ntype: nucleotide type, RNA or DNA
    :param double_stranded: is double stranded?
    :return: None
    """
    df["mw"] = df["sequence"].apply(
        lambda x: sequence.get_molecular_weight(x, ntype, double_stranded)
    )
    return df


def get_default_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds names to dataframe, if not already present
    :param df: dataframe
    :return: dataframe with names
    """
    if "name" in df.columns:
        raise ValueError("Dataframe already has names")
    # add `name` column to dataframe in the form of `seq_1`, `seq_2`, etc.
    df["name"] = df.index.map(lambda x: "seq_" + str(x))
    return df


def get_reverse_complement(df: pd.DataFrame, ntype: str) -> pd.DataFrame:
    """
    reverse complements each sequence in the dataframe
    :param df: dataframe
    :param ntype: nucleotide type, RNA or DNA
    :return: stores reverse complement in dataframe rev_comp column
    """
    df["rev_comp"] = df["sequence"].apply(
        lambda x: sequence.get_reverse_complement(x, ntype)
    )
    return df


def has_5p_sequence(df: pd.DataFrame, p5_seq: str) -> bool:
    """
    checks to see if p5_seq is present in the 5' end of the sequence
    :param df: dataframe
    :return: True if 5' sequence is present, False otherwise
    """
    return df["sequence"].str.startswith(p5_seq).all()


def has_3p_sequence(df: pd.DataFrame, p3_seq: str) -> bool:
    """
    checks to see if p5_seq is present in the 3' end of the sequence
    :param df: dataframe
    :return: True if 3' sequence is present, False otherwise
    """
    return df["sequence"].str.endswith(p3_seq).all()


def has_sequence(df: pd.DataFrame, seq: str) -> bool:
    """
    checks to see if seq is present in the sequence
    :param df: dataframe
    :return: True if sequence is present, False otherwise
    """
    return df["sequence"].str.contains(seq).all()


def has_t7_promoter(df: pd.DataFrame) -> bool:
    """
    checks if each sequence in the dataframe has a T7 promoter
    :param df: dataframe
    :return: None
    """
    has_t7 = df[df["sequence"].str.startswith("TTCTAATACGACTCACTATA")]
    if len(has_t7) != len(df):
        return False
    return True


def to_dna(df: pd.DataFrame) -> pd.DataFrame:
    """
    converts each sequence in dataframe to DNA
    :return: None
    """
    df["sequence"] = df["sequence"].apply(sequence.to_dna)
    if "structure" in df.columns:
        df = df.drop(columns=["structure"])
    return df


def to_dna_template(df: pd.DataFrame) -> pd.DataFrame:
    """
    converts each sequence in dataframe to DNA
    :return: None
    """
    df["sequence"] = df["sequence"].apply(sequence.to_dna_template)
    if "structure" in df.columns:
        df = df.drop(columns=["structure"])
    return df


def to_fasta(df: pd.DataFrame, filename: str) -> None:
    """
    writes the sequences in the dataframe to a fasta file
    :param df: dataframe
    :param filename: fasta file path
    :return: None
    """
    with open(filename, "w", encoding="utf-8") as f:
        for _, row in df.iterrows():
            f.write(f">{row['name']}\n")
            f.write(f"{row['sequence']}\n")


def to_opool(df: pd.DataFrame, name: str, filename: str) -> None:
    """
    writes the sequences in the dataframe to an opool file
    :param df: dataframe
    :param name: opool name
    :param filename: opool file path
    :return: None
    """
    df["name"] = name
    df = df[["name", "sequence"]]
    df.to_xlsx(filename, index=False)


def to_rna(df: pd.DataFrame) -> pd.DataFrame:
    """
    converts each sequence in dataframe to DNA
    :return: None
    """
    df["sequence"] = df["sequence"].apply(sequence.to_rna)
    return df


def trim(df, p5_length, p3_length) -> pd.DataFrame:
    """
    takes a data frame and trims the sequences. If there is a structure
    it will also trim the structure
    :param df: dataframe
    :param p5_length: length to trim from 5'
    :param p3_length: length to trim from 3'
    :return: None
    """
    # trim `sequence` column and `structure` column
    p3_length = -p3_length
    if p5_length == 0:
        p5_length = None
    if p3_length == 0:
        p3_length = None
    df["sequence"] = df["sequence"].str.slice(p5_length, p3_length)
    if "structure" in df.columns:
        df["structure"] = df["structure"].str.slice(p5_length, p3_length)
    return df


def transcribe(df: pd.DataFrame, ignore_missing_t7=False) -> pd.DataFrame:
    """
    transcribes each sequence in the dataframe (DNA -> RNA) removes t7 promoter
    :param df: dataframe with DNA template sequences
    :param ignore_missing_t7: ignore sequences that don't have a T7 promoter
    :return: dataframe with RNA sequences
    """
    if not has_t7_promoter(df) and not ignore_missing_t7:
        raise ValueError("not all sequences start with T7 promoter")
    if not ignore_missing_t7:
        df = trim(df, 20, 0)
    df = to_rna(df)
    df = fold(df)
    return df
