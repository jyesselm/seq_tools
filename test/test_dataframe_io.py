"""
Comprehensive tests for dataframe.io module
"""

import os

import pandas as pd
import pytest

from seq_tools.dataframe.io import to_opool


class TestToOpool:
    """Tests for to_opool function."""

    @pytest.mark.skipif(
        not hasattr(pd.DataFrame, "to_excel"), reason="openpyxl not installed"
    )
    def test_to_opool_basic(self):
        """Test basic opool export."""
        pytest.importorskip("openpyxl")
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["ATCGATCG", "GCTAGCTA"],
            }
        )
        filename = "test_opool.xlsx"
        try:
            to_opool(df, "my_pool", filename)
            assert os.path.exists(filename)
            # Read back and verify
            df_read = pd.read_excel(filename)
            assert len(df_read) == 2
            assert all(df_read["name"] == "my_pool")
        finally:
            if os.path.exists(filename):
                os.remove(filename)

    @pytest.mark.skipif(
        not hasattr(pd.DataFrame, "to_excel"), reason="openpyxl not installed"
    )
    def test_to_opool_with_extra_columns(self):
        """Test opool export with extra columns that should be removed."""
        pytest.importorskip("openpyxl")
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["ATCGATCG", "GCTAGCTA"],
                "extra": ["val1", "val2"],
            }
        )
        filename = "test_opool_extra.xlsx"
        try:
            to_opool(df, "my_pool", filename)
            assert os.path.exists(filename)
            df_read = pd.read_excel(filename)
            # Should only have name and sequence columns
            assert list(df_read.columns) == ["name", "sequence"]
        finally:
            if os.path.exists(filename):
                os.remove(filename)

    @pytest.mark.skipif(
        not hasattr(pd.DataFrame, "to_excel"), reason="openpyxl not installed"
    )
    def test_to_opool_single_sequence(self):
        """Test opool export with single sequence."""
        pytest.importorskip("openpyxl")
        df = pd.DataFrame(
            {
                "name": ["seq1"],
                "sequence": ["ATCGATCG"],
            }
        )
        filename = "test_opool_single.xlsx"
        try:
            to_opool(df, "my_pool", filename)
            assert os.path.exists(filename)
            df_read = pd.read_excel(filename)
            assert len(df_read) == 1
        finally:
            if os.path.exists(filename):
                os.remove(filename)
