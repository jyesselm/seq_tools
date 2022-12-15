"""
module to test dataframe.py
"""
import pandas as pd
from seq_tools.dataframe import to_dna, to_dna_template, to_rna, add


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


def test_add():
    """
    test add function
    """
    df = get_test_data_dna()
    df = add(df, "AAAA", "CCCC")
    assert df["sequence"][0] == "AAAAGGGGTTTTCCCCCCCC"


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
