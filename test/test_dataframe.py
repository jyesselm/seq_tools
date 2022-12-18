"""
module to test dataframe.py
"""
import pandas as pd
from seq_tools.dataframe import (
    add,
    get_extinction_coeff,
    get_molecular_weight,
    fold,
    to_dna,
    to_dna_template,
    to_rna,
    trim,
    transcribe
)


# generate test data ################################################################


def get_test_data_rna() -> pd.DataFrame:
    """
    get test data for dna
    :return: pd.DataFrame
    """
    return pd.DataFrame(
        [
            ["seq_0", "GGGGUUUUCCCC", "((((....))))"],
        ],
        columns=["name", "sequence", "structure"],
    )


def get_test_data_dna() -> pd.DataFrame:
    """
    get test data for dna
    :return: pd.DataFrame
    """
    return pd.DataFrame(
        [
            ["seq_0", "GGGGTTTTCCCC"],
        ],
        columns=["name", "sequence"],
    )


# tests  ##########################################################################


def test_add():
    """
    test add function
    """
    df = get_test_data_dna()
    df = add(df, "AAAA", "CCCC")
    assert df["sequence"][0] == "AAAAGGGGTTTTCCCCCCCC"


def test_get_extinction_coeff_dna():
    """
    test get_extinction_coeff function
    """
    df = pd.DataFrame([["seq_0", "ACGT"]], columns=["name", "sequence"])
    df = to_dna(df)
    df = get_extinction_coeff(df, "DNA", True)
    assert df["extinction_coeff"][0] == 66656


def test_get_extinction_coeff_rna():
    """
    test single stranded RNA extinction coefficient
    """
    df = pd.DataFrame([["seq_0", "AAAAAAAAUUUU"]], columns=["name", "sequence"])
    df = get_extinction_coeff(df, "RNA", False)
    assert df["extinction_coeff"][0] == 137100


def test_get_extinction_coeff_rna_w_struc():
    """
    test structured RNA extinction coefficient
    """
    df = pd.DataFrame(
        [["seq_0", "AAAAAAAAUUUU", "((((....))))"]],
        columns=["name", "sequence", "structure"],
    )
    df = get_extinction_coeff(df, "RNA", False)
    assert df["extinction_coeff"][0] == 113336


def test_get_molecular_weight_rna():
    """
    test get_molecular_weight function
    """
    df = pd.DataFrame([["seq_0", "AUG"]], columns=["name", "sequence"])
    df = get_molecular_weight(df, "RNA", False)
    assert df["mw"][0] == 1034.6


def test_fold():
    """
    test fold function
    """
    df = get_test_data_rna()
    df = df[["name", "sequence"]]
    df = fold(df)
    assert df["structure"][0] == "((((....))))"


def test_to_dna():
    """
    test to_dna function
    """
    df = get_test_data_rna()
    df = to_dna(df)
    assert df["sequence"][0] == "GGGGTTTTCCCC"


def test_to_dna_template():
    """
    test to_dna_template function
    """
    df = get_test_data_rna()
    df = to_dna_template(df)
    assert df["sequence"][0] == "TTCTAATACGACTCACTATAGGGGTTTTCCCC"


def test_to_rna():
    """
    test to_rna function
    """
    df = get_test_data_dna()
    df = to_rna(df)
    assert df["sequence"][0] == "GGGGUUUUCCCC"


def test_trim():
    """
    test trim function
    """
    df = get_test_data_dna()
    df = trim(df, 1, 2)
    assert df["sequence"][0] == "GGGTTTTCC"
    # okay with this behavior basically removes the entire string if the trim is too large
    df = get_test_data_dna()
    df = trim(df, 100, 100)
    assert df["sequence"][0] == ""
    df = get_test_data_dna()
    df.iloc[0]["sequence"] = "TTCTAATACGACTCACTATAGGGGTTTTCCCC"


def test_transcribe():
    """
    test transcribe function
    """
    df = get_test_data_dna()
    df.iloc[0]["sequence"] = "TTCTAATACGACTCACTATAGGGGTTTTCCCC"
    df = transcribe(df)
    assert df["sequence"][0] == "GGGGUUUUCCCC"
