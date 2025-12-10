"""
Tests for utils module
"""

from pathlib import Path

import pandas as pd
import pytest

from seq_tools.utils import (
    dataframe_to_sequences,
    get_resources_path,
    sequence_to_dataframe,
    sequences_to_dataframe,
)


class TestGetResourcesPath:
    """Tests for get_resources_path function."""

    def test_get_resources_path_exists(self):
        """Test that resources path exists and is a directory."""
        path = get_resources_path()
        assert isinstance(path, Path)
        assert path.exists()
        assert path.is_dir()

    def test_get_resources_path_contains_expected_files(self):
        """Test that resources directory contains expected CSV files."""
        path = get_resources_path()
        # Check for expected resource files
        assert (path / "p5_sequences.csv").exists()
        assert (path / "p3_sequences.csv").exists()

    def test_get_resources_path_returns_same_path(self):
        """Test that multiple calls return the same path."""
        path1 = get_resources_path()
        path2 = get_resources_path()
        assert path1 == path2


class TestSequenceToDataframe:
    """Tests for sequence_to_dataframe function."""

    def test_sequence_to_dataframe_basic(self):
        """Test basic conversion of sequence to dataframe."""
        df = sequence_to_dataframe("ATCG", name="test_seq")
        assert len(df) == 1
        assert "name" in df.columns
        assert "sequence" in df.columns
        assert df.iloc[0]["name"] == "test_seq"
        assert df.iloc[0]["sequence"] == "ATCG"

    def test_sequence_to_dataframe_default_name(self):
        """Test conversion with default name."""
        df = sequence_to_dataframe("ATCG")
        assert len(df) == 1
        assert df.iloc[0]["name"] == "sequence"
        assert df.iloc[0]["sequence"] == "ATCG"

    def test_sequence_to_dataframe_long_sequence(self):
        """Test conversion of long sequence."""
        long_seq = "ATCG" * 100
        df = sequence_to_dataframe(long_seq, name="long_seq")
        assert len(df) == 1
        assert df.iloc[0]["sequence"] == long_seq

    def test_sequence_to_dataframe_empty_string(self):
        """Test conversion of empty sequence string."""
        df = sequence_to_dataframe("", name="empty")
        assert len(df) == 1
        assert df.iloc[0]["sequence"] == ""

    def test_sequence_to_dataframe_rna(self):
        """Test conversion of RNA sequence."""
        df = sequence_to_dataframe("AUCG", name="rna_seq")
        assert df.iloc[0]["sequence"] == "AUCG"


class TestSequencesToDataframe:
    """Tests for sequences_to_dataframe function."""

    def test_sequences_to_dataframe_basic(self):
        """Test basic conversion of multiple sequences."""
        sequences = ["ATCG", "GCTA", "TTAA"]
        names = ["seq1", "seq2", "seq3"]
        df = sequences_to_dataframe(sequences, names=names)
        assert len(df) == 3
        assert list(df["name"]) == names
        assert list(df["sequence"]) == sequences

    def test_sequences_to_dataframe_default_names(self):
        """Test conversion with auto-generated names."""
        sequences = ["ATCG", "GCTA", "TTAA"]
        df = sequences_to_dataframe(sequences)
        assert len(df) == 3
        assert list(df["name"]) == ["seq_0", "seq_1", "seq_2"]
        assert list(df["sequence"]) == sequences

    def test_sequences_to_dataframe_single_sequence(self):
        """Test conversion of single sequence in list."""
        sequences = ["ATCG"]
        df = sequences_to_dataframe(sequences)
        assert len(df) == 1
        assert df.iloc[0]["name"] == "seq_0"
        assert df.iloc[0]["sequence"] == "ATCG"

    def test_sequences_to_dataframe_empty_list(self):
        """Test conversion of empty sequence list."""
        sequences = []
        df = sequences_to_dataframe(sequences)
        assert len(df) == 0
        assert "name" in df.columns
        assert "sequence" in df.columns

    def test_sequences_to_dataframe_mismatched_lengths(self):
        """Test that mismatched sequence and name counts raise error."""
        sequences = ["ATCG", "GCTA", "TTAA"]
        names = ["seq1", "seq2"]  # Too few names
        with pytest.raises(ValueError, match="must match"):
            sequences_to_dataframe(sequences, names=names)

    def test_sequences_to_dataframe_mismatched_lengths_too_many_names(self):
        """Test error with too many names."""
        sequences = ["ATCG", "GCTA"]
        names = ["seq1", "seq2", "seq3"]  # Too many names
        with pytest.raises(ValueError, match="must match"):
            sequences_to_dataframe(sequences, names=names)

    def test_sequences_to_dataframe_mixed_lengths(self):
        """Test conversion with sequences of different lengths."""
        sequences = ["ATCG", "GC", "TTAATTAA"]
        df = sequences_to_dataframe(sequences)
        assert len(df) == 3
        assert df.iloc[0]["sequence"] == "ATCG"
        assert df.iloc[1]["sequence"] == "GC"
        assert df.iloc[2]["sequence"] == "TTAATTAA"

    def test_sequences_to_dataframe_rna_sequences(self):
        """Test conversion of RNA sequences."""
        sequences = ["AUCG", "GCUA", "UUAA"]
        names = ["rna1", "rna2", "rna3"]
        df = sequences_to_dataframe(sequences, names=names)
        assert len(df) == 3
        assert list(df["sequence"]) == sequences


class TestDataframeToSequences:
    """Tests for dataframe_to_sequences function."""

    def test_dataframe_to_sequences_basic(self):
        """Test basic conversion of dataframe to sequence list."""
        df = pd.DataFrame({"name": ["seq1", "seq2"], "sequence": ["ATCG", "GCTA"]})
        result = dataframe_to_sequences(df)
        assert isinstance(result, list)
        assert len(result) == 2
        assert result[0] == {"name": "seq1", "sequence": "ATCG"}
        assert result[1] == {"name": "seq2", "sequence": "GCTA"}

    def test_dataframe_to_sequences_with_extra_columns(self):
        """Test conversion includes all columns."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["ATCG", "GCTA"],
                "structure": ["(..)", "(..)"],
                "mfe": [-1.0, -2.0],
            }
        )
        result = dataframe_to_sequences(df)
        assert len(result) == 2
        assert result[0] == {
            "name": "seq1",
            "sequence": "ATCG",
            "structure": "(..)",
            "mfe": -1.0,
        }

    def test_dataframe_to_sequences_single_row(self):
        """Test conversion of single row dataframe."""
        df = pd.DataFrame({"name": ["seq1"], "sequence": ["ATCG"]})
        result = dataframe_to_sequences(df)
        assert len(result) == 1
        assert result[0] == {"name": "seq1", "sequence": "ATCG"}

    def test_dataframe_to_sequences_empty(self):
        """Test conversion of empty dataframe."""
        df = pd.DataFrame({"name": [], "sequence": []})
        result = dataframe_to_sequences(df)
        assert isinstance(result, list)
        assert len(result) == 0

    def test_dataframe_to_sequences_preserves_types(self):
        """Test that data types are preserved in conversion."""
        df = pd.DataFrame(
            {
                "name": ["seq1"],
                "sequence": ["ATCG"],
                "length": [4],
                "mfe": [-1.5],
                "flag": [True],
            }
        )
        result = dataframe_to_sequences(df)
        assert result[0]["length"] == 4
        assert result[0]["mfe"] == -1.5
        assert result[0]["flag"] is True


class TestRoundTripConversion:
    """Tests for round-trip conversions between formats."""

    def test_sequence_to_dataframe_to_sequences(self):
        """Test round-trip: sequence -> dataframe -> sequences."""
        original_seq = "ATCGATCG"
        original_name = "test_seq"

        # sequence -> dataframe
        df = sequence_to_dataframe(original_seq, name=original_name)

        # dataframe -> sequences
        result = dataframe_to_sequences(df)

        assert len(result) == 1
        assert result[0]["sequence"] == original_seq
        assert result[0]["name"] == original_name

    def test_sequences_to_dataframe_to_sequences(self):
        """Test round-trip: sequences -> dataframe -> sequences."""
        original_sequences = ["ATCG", "GCTA", "TTAA"]
        original_names = ["seq1", "seq2", "seq3"]

        # sequences -> dataframe
        df = sequences_to_dataframe(original_sequences, names=original_names)

        # dataframe -> sequences
        result = dataframe_to_sequences(df)

        assert len(result) == 3
        for i, item in enumerate(result):
            assert item["sequence"] == original_sequences[i]
            assert item["name"] == original_names[i]
