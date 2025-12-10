"""
Comprehensive tests for dataframe.analysis module
"""

import pandas as pd

from seq_tools.dataframe.analysis import (
    calc_edit_distance,
    calc_edit_distance_parallel,
    determine_ntype,
)


class TestCalcEditDistance:
    """Tests for calc_edit_distance function."""

    def test_calc_edit_distance_identical_sequences(self):
        """Test edit distance with identical sequences."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2", "seq3"],
                "sequence": ["ATCGATCG", "ATCGATCG", "ATCGATCG"],
            }
        )
        result = calc_edit_distance(df)
        assert result == 0.0

    def test_calc_edit_distance_completely_different(self):
        """Test edit distance with completely different sequences."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["AAAAAAAA", "TTTTTTTT"],
            }
        )
        result = calc_edit_distance(df)
        assert result == 8.0

    def test_calc_edit_distance_multiple_sequences(self):
        """Test edit distance with multiple sequences."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2", "seq3", "seq4"],
                "sequence": ["ATCGATCG", "ATCGATCG", "ATCGATCC", "GGGGGGGG"],
            }
        )
        result = calc_edit_distance(df)
        assert isinstance(result, float)
        assert result >= 0


class TestCalcEditDistanceParallel:
    """Tests for calc_edit_distance_parallel function."""

    def test_calc_edit_distance_parallel_with_processes(self):
        """Test parallel edit distance with processes."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2", "seq3"],
                "sequence": ["ATCGATCG", "ATCGATCG", "GCTAGCTA"],
            }
        )
        result = calc_edit_distance_parallel(df, n_workers=2, use_threads=False)
        assert isinstance(result, float)
        assert result >= 0

    def test_calc_edit_distance_parallel_custom_workers(self):
        """Test parallel edit distance with custom number of workers."""
        df = pd.DataFrame(
            {
                "name": [f"seq_{i}" for i in range(10)],
                "sequence": ["ATCGATCG"] * 5 + ["GCTAGCTA"] * 5,
            }
        )
        result = calc_edit_distance_parallel(df, n_workers=4, use_threads=True)
        assert isinstance(result, float)
        assert result >= 0

    def test_calc_edit_distance_parallel_identical(self):
        """Test parallel edit distance with identical sequences."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["ATCGATCG", "ATCGATCG"],
            }
        )
        result = calc_edit_distance_parallel(df, use_threads=True)
        assert result == 0.0


class TestDetermineNtype:
    """Tests for determine_ntype function."""

    def test_determine_ntype_all_dna(self):
        """Test determining DNA type."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["ATCGATCGATCG", "TTTTAAAACCCC"],
            }
        )
        assert determine_ntype(df) == "DNA"

    def test_determine_ntype_all_rna(self):
        """Test determining RNA type."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["AUCGAUCGAUCG", "UUUUAAAACCCC"],
            }
        )
        assert determine_ntype(df) == "RNA"

    def test_determine_ntype_short_ambiguous(self):
        """Test determining type with short ambiguous sequences."""
        df = pd.DataFrame({"name": ["seq1", "seq2"], "sequence": ["AGCAGC", "CGACGA"]})
        # Should default to DNA for ambiguous cases
        result = determine_ntype(df)
        assert result in ["DNA", "RNA"]

    def test_determine_ntype_mixed_short(self):
        """Test with mixed T and U in short sequences."""
        df = pd.DataFrame({"name": ["seq1", "seq2"], "sequence": ["ATCG", "AUCG"]})
        # Short sequences should not raise error
        result = determine_ntype(df)
        assert result in ["DNA", "RNA"]

    def test_determine_ntype_single_sequence_dna(self):
        """Test with single DNA sequence."""
        df = pd.DataFrame({"name": ["seq1"], "sequence": ["ATCGATCGATCGATCG"]})
        assert determine_ntype(df) == "DNA"

    def test_determine_ntype_single_sequence_rna(self):
        """Test with single RNA sequence."""
        df = pd.DataFrame({"name": ["seq1"], "sequence": ["AUCGAUCGAUCGAUCG"]})
        assert determine_ntype(df) == "RNA"
