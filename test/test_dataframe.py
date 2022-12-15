"""
module to test dataframe.py
"""
import pandas as pd
from seq_tools.dataframe import to_dna


def get_test_data_rna() -> pd.DataFrame:
    """
    get test data for dna
    :return: pd.DataFrame
    """
    return pd.DataFrame(
        [
            ["seq_0", "AUGCCG"],
        ],
        columns=["name", "sequence"],
    )


def test_to_dna():
    """
    test to_dna function
    :return: None
    """
    df = get_test_data_rna()
    to_dna(df)
    assert df["sequence"][0] == "ATGCCG"
