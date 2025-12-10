"""
Additional comprehensive tests for dataframe.core module
"""

import pandas as pd

from seq_tools.dataframe.core import (
    fold,
    run_in_parallel,
    split,
    trim,
)


class TestTrimEdgeCases:
    """Additional tests for trim function edge cases."""

    def test_trim_with_list_column(self):
        """Test trimming columns containing lists."""
        df = pd.DataFrame(
            {
                "name": ["seq1"],
                "sequence": ["AAAATCGGGG"],
                "quality": [[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]],
            }
        )
        result = trim(df, 2, 2, extra_columns=["quality"])
        assert result["sequence"].iloc[0] == "AATCGG"
        assert result["quality"].iloc[0] == [30, 40, 50, 60, 70, 80]

    def test_trim_both_zero(self):
        """Test trim with both start and end as 0."""
        df = pd.DataFrame(
            {
                "name": ["seq1"],
                "sequence": ["ATCGATCG"],
            }
        )
        result = trim(df, 0, 0)
        assert result["sequence"].iloc[0] == "ATCGATCG"

    def test_trim_missing_extra_column(self):
        """Test trim when extra column doesn't exist."""
        df = pd.DataFrame(
            {
                "name": ["seq1"],
                "sequence": ["ATCGATCG"],
            }
        )
        # Should not raise error, just skip missing column
        result = trim(df, 2, 2, extra_columns=["nonexistent"])
        assert result["sequence"].iloc[0] == "CGAT"


class TestSplitEdgeCases:
    """Additional tests for split function edge cases."""

    def test_split_single_chunk(self):
        """Test split with single chunk."""
        df = pd.DataFrame(
            {"name": [f"seq_{i}" for i in range(5)], "sequence": ["ATCG"] * 5}
        )
        chunks = split(df, 1)
        assert len(chunks) == 1
        assert len(chunks[0]) == 5

    def test_split_more_chunks_than_rows(self):
        """Test split with more chunks than rows."""
        df = pd.DataFrame({"name": ["seq_0", "seq_1"], "sequence": ["ATCG", "GCTA"]})
        chunks = split(df, 5)
        # Should create empty chunks
        assert len(chunks) == 6


class TestFoldFunction:
    """Additional tests for fold function."""

    def test_fold_multiple_sequences(self):
        """Test folding multiple sequences."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["GGGGUUUUCCCC", "AAAAUUUU"],
            }
        )
        result = fold(df)
        assert "structure" in result.columns
        assert "mfe" in result.columns
        assert "ens_defect" in result.columns
        assert len(result) == 2


class TestRunInParallelEdgeCases:
    """Additional tests for run_in_parallel function."""

    def test_run_in_parallel_single_thread(self):
        """Test parallel execution with single thread."""
        df = pd.DataFrame(
            {"name": [f"seq_{i}" for i in range(5)], "sequence": ["ATCG"] * 5}
        )

        def add_column(chunk_df):
            chunk_df = chunk_df.copy()
            chunk_df["test"] = "value"
            return chunk_df

        result = run_in_parallel(df, add_column, 1)
        assert len(result) == 5
        assert "test" in result.columns
