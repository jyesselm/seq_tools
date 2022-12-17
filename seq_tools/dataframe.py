"""
module for working with dataframes that contain nucleotide sequences
"""
import pandas as pd
import vienna

from seq_tools import sequence, extinction_coeff


def add(df, p5_seq, p3_seq) -> pd.DataFrame:
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


def get_extinction_coeff(df, ntype, double_stranded):
    """
    calculates the extinction coefficient for each sequence in the dataframe
    :param df: dataframe
    :param ntype: nucleotide type, RNA or DNA
    :param double_stranded: is double stranded?
    :return: None
    """
    if ntype == "RNA" and "structure" in df.columns:
        df["extinction_coefficient"] = df.apply(
            lambda x: extinction_coeff.get_extinction_coeff(
                x["sequence"], ntype, double_stranded, x["structure"]
            )
        )
    else:
        df["extinction_coeff"] = df["sequence"].apply(
            lambda x: extinction_coeff.get_extinction_coeff(
                x, ntype, double_stranded
            )
        )
    return df


def get_molecular_weight(df, ntype, double_stranded) -> pd.DataFrame:
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


def fold(df) -> pd.DataFrame:
    """
    folds each sequence in the dataframe
    :param df: dataframe
    :return: None
    """

    def _fold(seq):
        v_res = vienna.fold(seq)
        return pd.Series(
            [v_res.dot_bracket, v_res.mfe, v_res.ens_defect],
            index=["structure", "mfe", "ens_defect"],
        )

    df[["structure", "mfe", "ens_defect"]] = df["sequence"].apply(_fold)
    return df


def to_dna(df) -> pd.DataFrame:
    """
    converts each sequence in dataframe to DNA
    :return: None
    """
    df["sequence"] = df["sequence"].apply(sequence.to_dna)
    if "structure" in df.columns:
        df = df.drop(columns=["structure"])
    return df


def to_dna_template(df) -> pd.DataFrame:
    """
    converts each sequence in dataframe to DNA
    :return: None
    """
    df["sequence"] = df["sequence"].apply(sequence.to_dna_template)
    if "structure" in df.columns:
        df = df.drop(columns=["structure"])
    return df


def to_rna(df) -> pd.DataFrame:
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
    df["sequence"] = df["sequence"].str.slice(p5_length, p3_length)
    if "structure" in df.columns:
        df["structure"] = df["structure"].str.slice(p5_length, p3_length)
    return df


def transcribe(df):
    """
    transcribes each sequence in the dataframe (DNA -> RNA)
    :param df:
    :return:
    """
    has_t7 = df[df["sequence"].str.startswith("TTCTAATACGACTCACTATA")]
    if len(has_t7) != len(df):
        raise ValueError("not all sequences start with T7 promoter")
    df = trim(df, 20, 0)
    df = fold(df)
    return df
