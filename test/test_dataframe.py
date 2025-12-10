"""
module to test dataframe.py
"""

import os

import pandas as pd
import pytest

from seq_tools.dataframe import (
    add,
    calc_edit_distance,
    determine_ntype,
    fold,
    get_extinction_coeff,
    get_molecular_weight,
    get_reverse_complement,
    has_3p_sequence,
    has_5p_sequence,
    has_t7_promoter,
    to_dna,
    to_dna_template,
    to_fasta,
    to_rna,
    transcribe,
    trim,
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
    df.loc[0, "sequence"] = "TTCTAATACGACTCACTATA" + "GAAAATTTTGGGGCCCC"
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
    with open("test.fasta", encoding="utf-8") as f:
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
    # okay with this behavior basically removes the entire string if the trim is
    # too large
    df = get_test_data_dna()
    df = trim(df, 100, 100)
    assert df["sequence"][0] == ""
    df = get_test_data_dna()
    df.loc[0, "sequence"] = "TTCTAATACGACTCACTATAGGGGTTTTCCCC"


def test_transcribe():
    """
    test transcribe function
    """
    df = get_test_data_dna()
    df.loc[0, "sequence"] = "TTCTAATACGACTCACTATAGGGGTTTTCCCC"
    df = transcribe(df)
    assert df["sequence"][0] == "GGGGUUUUCCCC"


# Additional tests for core.py functions ##########################################


def test_split():
    """Test split function from dataframe.core."""
    from seq_tools.dataframe.core import split

    df = pd.DataFrame(
        {
            "name": [f"seq_{i}" for i in range(10)],
            "sequence": ["ATCG"] * 10,
        }
    )
    chunks = split(df, 3)
    assert len(chunks) == 4  # 3 equal chunks + 1 remainder
    assert len(chunks[0]) == 3
    assert len(chunks[1]) == 3
    assert len(chunks[2]) == 3
    assert len(chunks[3]) == 1


def test_split_even_division():
    """Test split with evenly divisible number of rows."""
    from seq_tools.dataframe.core import split

    df = pd.DataFrame(
        {
            "name": [f"seq_{i}" for i in range(9)],
            "sequence": ["ATCG"] * 9,
        }
    )
    chunks = split(df, 3)
    assert len(chunks) == 3
    assert all(len(chunk) == 3 for chunk in chunks)


def test_run_in_parallel():
    """Test run_in_parallel function from dataframe.core."""
    from seq_tools.dataframe.core import run_in_parallel

    df = pd.DataFrame(
        {
            "name": [f"seq_{i}" for i in range(10)],
            "sequence": ["ATCG"] * 10,
        }
    )

    def add_length_column(chunk_df):
        chunk_df = chunk_df.copy()
        chunk_df["length"] = chunk_df["sequence"].apply(len)
        return chunk_df

    result = run_in_parallel(df, add_length_column, 2)
    assert len(result) == 10
    assert "length" in result.columns
    assert all(result["length"] == 4)


def test_get_length():
    """Test get_length function."""
    from seq_tools.dataframe.core import get_length

    df = pd.DataFrame(
        {
            "name": ["seq1", "seq2", "seq3"],
            "sequence": ["ATCG", "ATCGATCG", "AT"],
        }
    )
    result = get_length(df)
    assert "length" in result.columns
    assert list(result["length"]) == [4, 8, 2]


def test_get_default_names():
    """Test get_default_names function."""
    from seq_tools.dataframe.core import get_default_names

    df = pd.DataFrame({"sequence": ["ATCG", "GCTA", "TTAA"]})
    result = get_default_names(df)
    assert "name" in result.columns
    assert list(result["name"]) == ["seq_0", "seq_1", "seq_2"]


def test_get_default_names_raises_if_name_exists():
    """Test that get_default_names raises error if name column already exists."""
    from seq_tools.dataframe.core import get_default_names

    df = pd.DataFrame({"name": ["s1", "s2"], "sequence": ["ATCG", "GCTA"]})
    with pytest.raises(ValueError, match="already has names"):
        get_default_names(df)


def test_trim_with_structure():
    """Test trim function with structure column."""
    df = pd.DataFrame(
        {
            "name": ["seq1"],
            "sequence": ["AAAATCGGGGTTTT"],
            "structure": ["....(...)...."],
        }
    )
    result = trim(df, 4, 4)
    assert result["sequence"].iloc[0] == "TCGGGG"
    assert result["structure"].iloc[0] == "(...)"


def test_trim_no_end():
    """Test trim with end=0."""
    df = pd.DataFrame(
        {
            "name": ["seq1"],
            "sequence": ["AAAATCGGGG"],
        }
    )
    result = trim(df, 4, 0)
    assert result["sequence"].iloc[0] == "TCGGGG"


def test_trim_no_start():
    """Test trim with start=0."""
    df = pd.DataFrame(
        {
            "name": ["seq1"],
            "sequence": ["AAAATCGGGG"],
        }
    )
    result = trim(df, 0, 4)
    assert result["sequence"].iloc[0] == "AAAATC"


def test_trim_extra_columns():
    """Test trim with extra columns."""
    df = pd.DataFrame(
        {
            "name": ["seq1"],
            "sequence": ["AAAATCGGGG"],
            "custom": ["XXXXXXXXXX"],
        }
    )
    result = trim(df, 2, 2, extra_columns=["custom"])
    assert result["sequence"].iloc[0] == "AATCGG"
    assert result["custom"].iloc[0] == "XXXXXX"


def test_add_with_structure():
    """Test add function refolds when structure column exists."""
    df = pd.DataFrame(
        {
            "name": ["seq1"],
            "sequence": ["GGGGUUUUCCCC"],
            "structure": ["((((....))))"],
        }
    )
    result = add(df, "AA", "UU")
    assert result["sequence"].iloc[0] == "AAGGGGUUUUCCCCUU"
    # Structure should be refolded
    assert "structure" in result.columns
    assert "mfe" in result.columns


# Additional tests for analysis.py functions ######################################


def test_calc_edit_distance_parallel():
    """Test parallel edit distance calculation."""
    from seq_tools.dataframe.analysis import calc_edit_distance_parallel

    df = pd.DataFrame(
        {
            "name": ["seq1", "seq2", "seq3"],
            "sequence": ["ATCGATCG", "ATCGATCG", "GCTAGCTA"],
        }
    )
    # Two sequences are identical, one is different
    result = calc_edit_distance_parallel(df, n_workers=2, use_threads=True)
    assert isinstance(result, float)
    assert result >= 0


def test_calc_edit_distance_parallel_single_sequence():
    """Test parallel edit distance with single sequence."""
    from seq_tools.dataframe.analysis import calc_edit_distance_parallel

    df = pd.DataFrame({"name": ["seq1"], "sequence": ["ATCG"]})
    result = calc_edit_distance_parallel(df)
    assert result == 0.0


def test_determine_ntype_uncertain():
    """Test determine_ntype with short sequences that are uncertain."""
    df = pd.DataFrame(
        {
            "name": ["seq1", "seq2"],
            "sequence": ["AGCN", "AGCN"],  # Short, no T or U
        }
    )
    result = determine_ntype(df)
    # Should default to DNA for uncertain short sequences
    assert result == "DNA"


def test_determine_ntype_mixed_long_sequences():
    """Test determine_ntype raises error with mixed long sequences."""
    df = pd.DataFrame(
        {
            "name": ["seq1", "seq2"],
            "sequence": [
                "ATCGATCGATCGATCG",  # Long DNA
                "AUCGAUCGAUCGAUCG",  # Long RNA
            ],
        }
    )
    with pytest.raises(ValueError, match="Cannot determine nucleotide type"):
        determine_ntype(df)
