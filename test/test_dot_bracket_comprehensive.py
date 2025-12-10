"""
Comprehensive tests for dot_bracket module to improve coverage
"""

import pytest

from seq_tools.dot_bracket import dotbracket_to_pairtable, inverse_brackets


class TestInverseBrackets:
    """Tests for inverse_brackets function."""

    def test_inverse_brackets_basic(self):
        """Test basic inverse bracket mapping."""
        result = inverse_brackets("([{")
        assert result["("] == 0
        assert result["["] == 1
        assert result["{"] == 2

    def test_inverse_brackets_right(self):
        """Test inverse bracket mapping for right brackets."""
        result = inverse_brackets(")]}")
        assert result[")"] == 0
        assert result["]"] == 1
        assert result["}"] == 2


class TestDotbracketToPairtable:
    """Tests for dotbracket_to_pairtable function."""

    def test_dotbracket_to_pairtable_empty_raises_error(self):
        """Test that empty structure raises error."""
        with pytest.raises(ValueError, match="Cannot convert empty structure"):
            dotbracket_to_pairtable("")

    def test_dotbracket_to_pairtable_too_many_closing_brackets(self):
        """Test error with unbalanced closing brackets."""
        with pytest.raises(ValueError, match="Too many closing brackets"):
            dotbracket_to_pairtable("(())))")

    def test_dotbracket_to_pairtable_too_many_opening_brackets(self):
        """Test error with unbalanced opening brackets."""
        with pytest.raises(ValueError, match="Too many opening brackets"):
            dotbracket_to_pairtable("(((()")

    def test_dotbracket_to_pairtable_multistrand(self):
        """Test conversion with multistrand structure."""
        struct = "(((((..&..)))))"
        pt = dotbracket_to_pairtable(struct)
        # Verify pair table structure
        assert isinstance(pt, list)
        assert len(pt) == len(struct) - struct.count("&")

    def test_dotbracket_to_pairtable_different_bracket_types(self):
        """Test conversion with different bracket types."""
        struct = "([{<>}])"
        pt = dotbracket_to_pairtable(struct)
        assert isinstance(pt, list)
        # First and last should be paired
        assert pt[0] == len(pt) - 1
        assert pt[-1] == 0

    def test_dotbracket_to_pairtable_nested_brackets(self):
        """Test conversion with nested brackets."""
        struct = "((([[[]]])))"
        pt = dotbracket_to_pairtable(struct)
        assert isinstance(pt, list)
        assert len(pt) == len(struct)

    def test_dotbracket_to_pairtable_dots_only(self):
        """Test conversion with only unpaired bases."""
        struct = "....."
        pt = dotbracket_to_pairtable(struct)
        # All should be -1 (unpaired)
        assert all(x == -1 for x in pt)
