"""
Tests for validation module
"""

import pandas as pd
import pytest

from seq_tools.validation import (
    ensure_name_column,
    validate_dataframe,
    validate_sequence,
)


class TestValidateSequence:
    """Tests for validate_sequence function."""

    def test_validate_sequence_valid_dna(self):
        """Test validation of valid DNA sequence."""
        validate_sequence("ATCG")
        validate_sequence("ATCGATCG")
        validate_sequence("atcg")  # lowercase should work

    def test_validate_sequence_valid_dna_with_n(self):
        """Test validation of DNA sequence with ambiguous N."""
        validate_sequence("ATCGN", ntype="DNA", allow_ambiguous=True)
        validate_sequence("NNNNATCG", ntype="DNA", allow_ambiguous=True)

    def test_validate_sequence_valid_dna_no_ambiguous(self):
        """Test validation of DNA sequence without ambiguous nucleotides."""
        validate_sequence("ATCG", ntype="DNA", allow_ambiguous=False)

    def test_validate_sequence_invalid_dna_with_n_not_allowed(self):
        """Test that N is invalid when ambiguous nucleotides not allowed."""
        with pytest.raises(ValueError, match="Invalid DNA sequence"):
            validate_sequence("ATCGN", ntype="DNA", allow_ambiguous=False)

    def test_validate_sequence_valid_rna(self):
        """Test validation of valid RNA sequence."""
        validate_sequence("AUCG", ntype="RNA")
        validate_sequence("AUCGAUCG", ntype="RNA")
        validate_sequence("aucg", ntype="RNA")  # lowercase should work

    def test_validate_sequence_valid_rna_with_n(self):
        """Test validation of RNA sequence with ambiguous N."""
        validate_sequence("AUCGN", ntype="RNA", allow_ambiguous=True)
        validate_sequence("NNNNAUCG", ntype="RNA", allow_ambiguous=True)

    def test_validate_sequence_valid_rna_no_ambiguous(self):
        """Test validation of RNA sequence without ambiguous nucleotides."""
        validate_sequence("AUCG", ntype="RNA", allow_ambiguous=False)

    def test_validate_sequence_invalid_rna_with_n_not_allowed(self):
        """Test that N is invalid when ambiguous nucleotides not allowed."""
        with pytest.raises(ValueError, match="Invalid RNA sequence"):
            validate_sequence("AUCGN", ntype="RNA", allow_ambiguous=False)

    def test_validate_sequence_invalid_dna_characters(self):
        """Test validation fails with invalid DNA characters."""
        with pytest.raises(ValueError, match="Invalid DNA sequence"):
            validate_sequence("ATCGU", ntype="DNA")
        with pytest.raises(ValueError, match="Invalid DNA sequence"):
            validate_sequence("ATCG123", ntype="DNA")
        with pytest.raises(ValueError, match="Invalid DNA sequence"):
            validate_sequence("ATCG-XYZ", ntype="DNA")

    def test_validate_sequence_invalid_rna_characters(self):
        """Test validation fails with invalid RNA characters."""
        with pytest.raises(ValueError, match="Invalid RNA sequence"):
            validate_sequence("AUCGT", ntype="RNA")
        with pytest.raises(ValueError, match="Invalid RNA sequence"):
            validate_sequence("AUCG123", ntype="RNA")
        with pytest.raises(ValueError, match="Invalid RNA sequence"):
            validate_sequence("AUCG-XYZ", ntype="RNA")

    def test_validate_sequence_empty(self):
        """Test validation fails with empty sequence."""
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            validate_sequence("")

    def test_validate_sequence_auto_detect_dna(self):
        """Test auto-detection of DNA sequence type."""
        validate_sequence("ATCG")  # Should auto-detect as DNA
        validate_sequence("ATCGATCGATCG")

    def test_validate_sequence_auto_detect_rna(self):
        """Test auto-detection of RNA sequence type."""
        validate_sequence("AUCG")  # Should auto-detect as RNA
        validate_sequence("AUCGAUCGAUCG")

    def test_validate_sequence_auto_detect_invalid(self):
        """Test auto-detection fails with invalid characters."""
        with pytest.raises(ValueError, match="Invalid sequence"):
            validate_sequence("ATCGXYZ")
        with pytest.raises(ValueError, match="Invalid sequence"):
            validate_sequence("123456")

    def test_validate_sequence_auto_detect_with_n(self):
        """Test auto-detection with ambiguous nucleotides."""
        validate_sequence("ATCGN", allow_ambiguous=True)
        validate_sequence("AUCGN", allow_ambiguous=True)


class TestValidateDataframe:
    """Tests for validate_dataframe function."""

    def test_validate_dataframe_valid(self):
        """Test validation of valid dataframe."""
        df = pd.DataFrame({"name": ["seq1", "seq2"], "sequence": ["ATCG", "GCTA"]})
        validate_dataframe(df)

    def test_validate_dataframe_valid_no_require_name(self):
        """Test validation without requiring name column."""
        df = pd.DataFrame({"sequence": ["ATCG", "GCTA"]})
        validate_dataframe(df, require_name=False)

    def test_validate_dataframe_missing_sequence(self):
        """Test validation fails when sequence column is missing."""
        df = pd.DataFrame({"name": ["seq1", "seq2"]})
        with pytest.raises(ValueError, match="must have a 'sequence' column"):
            validate_dataframe(df)

    def test_validate_dataframe_missing_name(self):
        """Test validation fails when name column is missing and required."""
        df = pd.DataFrame({"sequence": ["ATCG", "GCTA"]})
        with pytest.raises(ValueError, match="must have a 'name' column"):
            validate_dataframe(df, require_name=True)

    def test_validate_dataframe_with_extra_columns(self):
        """Test validation passes with extra columns."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["ATCG", "GCTA"],
                "structure": ["(..)", "(..)"],
                "mfe": [-1.0, -2.0],
            }
        )
        validate_dataframe(df)

    def test_validate_dataframe_empty(self):
        """Test validation of empty dataframe with correct columns."""
        df = pd.DataFrame({"name": [], "sequence": []})
        validate_dataframe(df)


class TestEnsureNameColumn:
    """Tests for ensure_name_column function."""

    def test_ensure_name_column_already_exists(self):
        """Test that existing name column is preserved."""
        df = pd.DataFrame({"name": ["seq1", "seq2"], "sequence": ["ATCG", "GCTA"]})
        result = ensure_name_column(df)
        assert "name" in result.columns
        assert list(result["name"]) == ["seq1", "seq2"]

    def test_ensure_name_column_adds_default_names(self):
        """Test that default names are added when missing."""
        df = pd.DataFrame({"sequence": ["ATCG", "GCTA", "TTAA"]})
        result = ensure_name_column(df)
        assert "name" in result.columns
        assert list(result["name"]) == ["seq_0", "seq_1", "seq_2"]

    def test_ensure_name_column_empty_dataframe(self):
        """Test handling of empty dataframe."""
        df = pd.DataFrame({"sequence": []})
        result = ensure_name_column(df)
        assert "name" in result.columns
        assert len(result["name"]) == 0

    def test_ensure_name_column_single_row(self):
        """Test handling of single row dataframe."""
        df = pd.DataFrame({"sequence": ["ATCG"]})
        result = ensure_name_column(df)
        assert "name" in result.columns
        assert list(result["name"]) == ["seq_0"]

    def test_ensure_name_column_preserves_other_columns(self):
        """Test that other columns are preserved."""
        df = pd.DataFrame({"sequence": ["ATCG", "GCTA"], "structure": ["(..)", "(..)"]})
        result = ensure_name_column(df)
        assert "sequence" in result.columns
        assert "structure" in result.columns
        assert list(result["sequence"]) == ["ATCG", "GCTA"]

    def test_ensure_name_column_returns_copy(self):
        """Test that function returns a copy and doesn't modify original."""
        df = pd.DataFrame({"sequence": ["ATCG", "GCTA"]})
        result = ensure_name_column(df)
        # Original should not have name column
        assert "name" not in df.columns
        # Result should have name column
        assert "name" in result.columns
