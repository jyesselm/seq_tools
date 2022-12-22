"""
module to test dataframe.py
"""
import os
import pytest
import pandas as pd
from seq_tools.dataframe import (
    add,
    calc_edit_distance,
    determine_ntype,
    fold,
    has_t7_promoter,
    has_5p_sequence,
    has_3p_sequence,
    get_extinction_coeff,
    get_molecular_weight,
    get_reverse_complement,
    to_dna,
    to_dna_template,
    to_fasta,
    to_rna,
    trim,
    transcribe,
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


def test_calc_edit_distance():
    """
    test calc_edit_distance function
    """
    # should return zero for only 1 sequence
    df = get_test_data_dna()
    val = calc_edit_distance(df)
    assert val == 0
    # should return 1 for 2 sequences
    df = get_test_data_dna()
    df.loc[1] = ["seq_1", "GGGGTATTCCCC"]
    val = calc_edit_distance(df)
    assert val == 1


def test_determine_ntype():
    """
    test determine_ntype function
    """
    df = get_test_data_dna()
    assert determine_ntype(df) == "DNA"
    df = get_test_data_rna()
    assert determine_ntype(df) == "RNA"
    df = pd.concat([get_test_data_dna(), get_test_data_rna()])
    with pytest.raises(ValueError):
        determine_ntype(df)


def test_fold():
    """
    test fold function
    """
    df = get_test_data_rna()
    df = df[["name", "sequence"]]
    df = fold(df)
    assert df["structure"][0] == "((((....))))"


def test_has_t7_promoter():
    """
    test has_t7_promoter function
    :return:
    """
    df = get_test_data_dna()
    has_t7 = has_t7_promoter(df)
    assert not has_t7
    df.iloc[0]["sequence"] = "TTCTAATACGACTCACTATA" + "GAAAATTTTGGGGCCCC"
    has_t7 = has_t7_promoter(df)
    assert has_t7


def test_has_5p_sequence():
    """
    test has_5p_sequence function
    """
    df = get_test_data_dna()
    df.loc[1] = ["seq_1", "GGGTTTTCCCC"]
    has_5p = has_5p_sequence(df, "GGG")
    assert has_5p
    has_5p = has_5p_sequence(df, "GGGG")
    assert not has_5p


def test_has_3p_sequence():
    """
    test has_3p_sequence function
    """
    df = get_test_data_dna()
    df.loc[1] = ["seq_1", "GGGTTTTCCCC"]
    has_3p = has_3p_sequence(df, "CCCC")
    assert has_3p
    has_3p = has_3p_sequence(df, "GGGG")
    assert not has_3p


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


def test_reverse_complement():
    """
    test reverse_complement function
    """
    df = get_test_data_dna()
    df = get_reverse_complement(df, "DNA")
    assert df["sequence"][0] == "GGGGTTTTCCCC"


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


def test_to_fasta():
    """
    test to_fasta function
    """
    df = get_test_data_dna()
    to_fasta(df, "test.fasta")
    assert os.path.isfile("test.fasta")
    with open("test.fasta", "r", encoding="utf-8") as f:
        lines = f.readlines()
    assert lines[0] == ">seq_0\n"
    assert lines[1] == "GGGGTTTTCCCC\n"
    os.remove("test.fasta")


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
