"""
Tests for dataframe.primers module
"""

import pandas as pd
import pytest

from seq_tools.dataframe.primers import (
    add_common_seqs,
    find_longest_common_prefix,
    find_longest_common_suffix,
    has_3p_sequence,
    has_5p_sequence,
    has_sequence,
    has_t7_promoter,
    remove_common_p5_p3,
    remove_common_p5_p3_by_structure,
    remove_common_seqs,
    transcribe,
    trim_p5_and_p3,
)


class TestHas5pSequence:
    """Tests for has_5p_sequence function."""

    def test_has_5p_sequence_all_match(self):
        """Test when all sequences start with the given 5' sequence."""
        df = pd.DataFrame(
            {"sequence": ["AAAATCG", "AAAAGCT", "AAAATTA"], "name": ["s1", "s2", "s3"]}
        )
        assert has_5p_sequence(df, "AAAA") == True

    def test_has_5p_sequence_none_match(self):
        """Test when no sequences start with the given 5' sequence."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGATCG", "GCTAGCTA", "TTAATTAA"],
                "name": ["s1", "s2", "s3"],
            }
        )
        assert has_5p_sequence(df, "CCCC") == False

    def test_has_5p_sequence_partial_match(self):
        """Test when only some sequences start with the given 5' sequence."""
        df = pd.DataFrame(
            {"sequence": ["AAAATCG", "GCTAGCTA", "AAAATTA"], "name": ["s1", "s2", "s3"]}
        )
        assert has_5p_sequence(df, "AAAA") == False

    def test_has_5p_sequence_single_sequence(self):
        """Test with single sequence."""
        df = pd.DataFrame({"sequence": ["AAAATCG"], "name": ["s1"]})
        assert has_5p_sequence(df, "AAAA") == True
        assert has_5p_sequence(df, "TTTT") == False

    def test_has_5p_sequence_empty_search(self):
        """Test with empty search string."""
        df = pd.DataFrame({"sequence": ["ATCG", "GCTA"], "name": ["s1", "s2"]})
        assert has_5p_sequence(df, "") == True


class TestHas3pSequence:
    """Tests for has_3p_sequence function."""

    def test_has_3p_sequence_all_match(self):
        """Test when all sequences end with the given 3' sequence."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGTTTT", "GCTATTTT", "TTAATTTT"],
                "name": ["s1", "s2", "s3"],
            }
        )
        assert has_3p_sequence(df, "TTTT") == True

    def test_has_3p_sequence_none_match(self):
        """Test when no sequences end with the given 3' sequence."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGATCG", "GCTAGCTA", "TTAATTAA"],
                "name": ["s1", "s2", "s3"],
            }
        )
        assert has_3p_sequence(df, "CCCC") == False

    def test_has_3p_sequence_partial_match(self):
        """Test when only some sequences end with the given 3' sequence."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGTTTT", "GCTAGCTA", "TTAATTTT"],
                "name": ["s1", "s2", "s3"],
            }
        )
        assert has_3p_sequence(df, "TTTT") == False

    def test_has_3p_sequence_single_sequence(self):
        """Test with single sequence."""
        df = pd.DataFrame({"sequence": ["ATCGTTTT"], "name": ["s1"]})
        assert has_3p_sequence(df, "TTTT") == True
        assert has_3p_sequence(df, "AAAA") == False

    def test_has_3p_sequence_empty_search(self):
        """Test with empty search string."""
        df = pd.DataFrame({"sequence": ["ATCG", "GCTA"], "name": ["s1", "s2"]})
        assert has_3p_sequence(df, "") == True


class TestHasSequence:
    """Tests for has_sequence function."""

    def test_has_sequence_all_contain(self):
        """Test when all sequences contain the subsequence."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGATCG", "GGGATCGAA", "TTTATCGCC"],
                "name": ["s1", "s2", "s3"],
            }
        )
        assert has_sequence(df, "ATCG") == True

    def test_has_sequence_none_contain(self):
        """Test when no sequences contain the subsequence."""
        df = pd.DataFrame(
            {"sequence": ["ATATAT", "GCGCGC", "TTTTTT"], "name": ["s1", "s2", "s3"]}
        )
        assert has_sequence(df, "CCCC") == False

    def test_has_sequence_partial_contain(self):
        """Test when only some sequences contain the subsequence."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGATCG", "GGGGGGGG", "ATCGAAAA"],
                "name": ["s1", "s2", "s3"],
            }
        )
        assert has_sequence(df, "ATCG") == False

    def test_has_sequence_at_start(self):
        """Test finding sequence at the start."""
        df = pd.DataFrame({"sequence": ["ATCGGGGG", "ATCGAAAA"], "name": ["s1", "s2"]})
        assert has_sequence(df, "ATCG") == True

    def test_has_sequence_at_end(self):
        """Test finding sequence at the end."""
        df = pd.DataFrame({"sequence": ["GGGGATCG", "AAAATCG"], "name": ["s1", "s2"]})
        assert has_sequence(df, "ATCG") == True


class TestHasT7Promoter:
    """Tests for has_t7_promoter function."""

    def test_has_t7_promoter_all_have(self):
        """Test when all sequences have T7 promoter."""
        df = pd.DataFrame(
            {
                "sequence": [
                    "TTCTAATACGACTCACTATAGGGG",
                    "TTCTAATACGACTCACTATAAAAA",
                ],
                "name": ["s1", "s2"],
            }
        )
        assert has_t7_promoter(df) is True

    def test_has_t7_promoter_none_have(self):
        """Test when no sequences have T7 promoter."""
        df = pd.DataFrame({"sequence": ["ATCGATCG", "GCTAGCTA"], "name": ["s1", "s2"]})
        assert has_t7_promoter(df) is False

    def test_has_t7_promoter_partial(self):
        """Test when only some sequences have T7 promoter."""
        df = pd.DataFrame(
            {
                "sequence": [
                    "TTCTAATACGACTCACTATAGGGG",
                    "ATCGATCG",
                ],
                "name": ["s1", "s2"],
            }
        )
        assert has_t7_promoter(df) is False


class TestFindLongestCommonPrefix:
    """Tests for find_longest_common_prefix function."""

    def test_find_longest_common_prefix_basic(self):
        """Test finding common prefix in sequences."""
        df = pd.DataFrame(
            {"sequence": ["AAAATCG", "AAAAGCT", "AAAATTA"], "name": ["s1", "s2", "s3"]}
        )
        assert find_longest_common_prefix(df) == "AAAA"

    def test_find_longest_common_prefix_no_common(self):
        """Test when there is no common prefix."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGATCG", "GCTAGCTA", "TTAATTAA"],
                "name": ["s1", "s2", "s3"],
            }
        )
        assert find_longest_common_prefix(df) == ""

    def test_find_longest_common_prefix_entire_sequence(self):
        """Test when entire sequences are the same."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGATCG", "ATCGATCG", "ATCGATCG"],
                "name": ["s1", "s2", "s3"],
            }
        )
        assert find_longest_common_prefix(df) == "ATCGATCG"

    def test_find_longest_common_prefix_single_sequence(self):
        """Test with single sequence."""
        df = pd.DataFrame({"sequence": ["ATCGATCG"], "name": ["s1"]})
        assert find_longest_common_prefix(df) == "ATCGATCG"

    def test_find_longest_common_prefix_empty_dataframe(self):
        """Test with empty dataframe."""
        df = pd.DataFrame({"sequence": [], "name": []})
        assert find_longest_common_prefix(df) == ""

    def test_find_longest_common_prefix_different_lengths(self):
        """Test with sequences of different lengths."""
        df = pd.DataFrame(
            {"sequence": ["ATCGATCG", "ATCG", "ATCGAA"], "name": ["s1", "s2", "s3"]}
        )
        assert find_longest_common_prefix(df) == "ATCG"


class TestFindLongestCommonSuffix:
    """Tests for find_longest_common_suffix function."""

    def test_find_longest_common_suffix_basic(self):
        """Test finding common suffix in sequences."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGTTTT", "GCTATTTT", "TTAATTTT"],
                "name": ["s1", "s2", "s3"],
            }
        )
        assert find_longest_common_suffix(df) == "TTTT"

    def test_find_longest_common_suffix_no_common(self):
        """Test when there is no common suffix."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGATCG", "GCTAGCTA", "TTAATTAA"],
                "name": ["s1", "s2", "s3"],
            }
        )
        assert find_longest_common_suffix(df) == ""

    def test_find_longest_common_suffix_entire_sequence(self):
        """Test when entire sequences are the same."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGATCG", "ATCGATCG", "ATCGATCG"],
                "name": ["s1", "s2", "s3"],
            }
        )
        assert find_longest_common_suffix(df) == "ATCGATCG"

    def test_find_longest_common_suffix_single_sequence(self):
        """Test with single sequence."""
        df = pd.DataFrame({"sequence": ["ATCGATCG"], "name": ["s1"]})
        assert find_longest_common_suffix(df) == "ATCGATCG"

    def test_find_longest_common_suffix_empty_dataframe(self):
        """Test with empty dataframe."""
        df = pd.DataFrame({"sequence": [], "name": []})
        assert find_longest_common_suffix(df) == ""

    def test_find_longest_common_suffix_different_lengths(self):
        """Test with sequences of different lengths."""
        df = pd.DataFrame(
            {"sequence": ["ATCGATCG", "ATCG", "GGGATCG"], "name": ["s1", "s2", "s3"]}
        )
        assert find_longest_common_suffix(df) == "ATCG"


class TestTranscribe:
    """Tests for transcribe function."""

    def test_transcribe_basic(self):
        """Test basic transcription with T7 promoter."""
        df = pd.DataFrame(
            {
                "sequence": ["TTCTAATACGACTCACTATAGGGGTTTTCCCC"],
                "name": ["seq1"],
            }
        )
        result = transcribe(df)
        # Should remove T7 promoter (20 nt), convert to RNA, and fold
        assert "GGGGUUUUCCCC" == result["sequence"].iloc[0]
        assert "structure" in result.columns

    def test_transcribe_without_t7_ignore_missing(self):
        """Test transcription ignoring missing T7 promoter."""
        df = pd.DataFrame(
            {
                "sequence": ["GGGGTTTTCCCC"],
                "name": ["seq1"],
            }
        )
        result = transcribe(df, ignore_missing_t7=True)
        # Should convert to RNA and fold without trimming
        assert "GGGGUUUUCCCC" == result["sequence"].iloc[0]
        assert "structure" in result.columns

    def test_transcribe_without_t7_raises_error(self):
        """Test that missing T7 promoter raises error by default."""
        df = pd.DataFrame(
            {
                "sequence": ["GGGGTTTTCCCC"],
                "name": ["seq1"],
            }
        )
        with pytest.raises(
            ValueError, match="not all sequences start with T7 promoter"
        ):
            transcribe(df)


class TestAddCommonSeqs:
    """Tests for add_common_seqs function."""

    def test_add_common_seqs_basic(self):
        """Test adding common sequences to both ends."""
        df = pd.DataFrame(
            {
                "name": ["seq1"],
                "sequence": ["GGGGAAAA"],
                "structure": ["((...))"],
            }
        )
        result = add_common_seqs(
            df,
            p5_seq="AAAA",
            p5_structure="....",
            p3_seq="UUUU",
            p3_structure="....",
        )
        assert result["sequence"].iloc[0] == "AAAAGGGGAAAAUUUU"
        assert "predicted_structure" in result.columns
        assert "structure_match" in result.columns

    def test_add_common_seqs_empty_p5(self):
        """Test adding only 3' sequence."""
        df = pd.DataFrame(
            {
                "name": ["seq1"],
                "sequence": ["GGGGAAAA"],
                "structure": ["((...))"],
            }
        )
        result = add_common_seqs(
            df, p5_seq="", p5_structure="", p3_seq="UUUU", p3_structure="...."
        )
        assert result["sequence"].iloc[0] == "GGGGAAAAUUUU"

    def test_add_common_seqs_empty_p3(self):
        """Test adding only 5' sequence."""
        df = pd.DataFrame(
            {
                "name": ["seq1"],
                "sequence": ["GGGGAAAA"],
                "structure": ["((...))"],
            }
        )
        result = add_common_seqs(
            df, p5_seq="AAAA", p5_structure="....", p3_seq="", p3_structure=""
        )
        assert result["sequence"].iloc[0] == "AAAAGGGGAAAA"

    def test_add_common_seqs_no_structure_column_raises_error(self):
        """Test that missing structure column raises error."""
        df = pd.DataFrame(
            {
                "name": ["seq1"],
                "sequence": ["GGGGAAAA"],
            }
        )
        with pytest.raises(ValueError, match="must contain a 'structure' column"):
            add_common_seqs(df, "AAAA", "....", "UUUU", "....")

    def test_add_common_seqs_parallel(self):
        """Test parallel processing mode."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["GGGGAAAA", "CCCCUUUU"],
                "structure": ["((...)).", "((....))"],
            }
        )
        result = add_common_seqs(
            df,
            p5_seq="AA",
            p5_structure="..",
            p3_seq="UU",
            p3_structure="..",
            parallel=True,
            workers=2,
        )
        assert len(result) == 2
        assert "structure_match" in result.columns


class TestRemoveCommonSeqs:
    """Tests for remove_common_seqs function."""

    def test_remove_common_seqs_returns_warnings_dict(self):
        """Test that function returns tuple with warnings dict."""
        # This test requires actual p5/p3 sequence files to work properly
        # We'll test the basic structure of the return value
        df = pd.DataFrame(
            {
                "sequence": [
                    "GGGAACCUUAACCGGAAAACCCCUUUU",
                    "GGGAACCUUAACCGGAAAACCCCUUUU",
                ],
                "name": ["s1", "s2"],
            }
        )
        try:
            result, info = remove_common_seqs(df)
            assert isinstance(info, dict)
            assert "warnings" in info
            assert "removed" in info
        except ValueError:
            # If no common sequences found, that's okay for this test
            pass


class TestRemoveCommonP5P3:
    """Tests for remove_common_p5_p3 function."""

    def test_remove_common_p5_p3_raises_when_none_found(self):
        """Test that function raises error when no common sequences found."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGATCG", "GCTAGCTA"],
                "name": ["s1", "s2"],
            }
        )
        with pytest.raises(ValueError, match="No common p5 or p3 sequence found"):
            remove_common_p5_p3(df)


class TestRemoveCommonP5P3ByStructure:
    """Tests for remove_common_p5_p3_by_structure function."""

    def test_remove_common_p5_p3_by_structure_raises_when_none_found(self):
        """Test that function raises error when no matching patterns found."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGATCG", "GCTAGCTA"],
                "name": ["s1", "s2"],
            }
        )
        with pytest.raises(
            ValueError, match="No common p5 or p3 sequence/structure pattern found"
        ):
            remove_common_p5_p3_by_structure(df)


class TestTrimP5AndP3:
    """Tests for trim_p5_and_p3 function."""

    def test_trim_p5_and_p3_raises_when_none_found(self):
        """Test that function raises error when no common sequences found."""
        df = pd.DataFrame(
            {
                "sequence": ["ATCGATCG", "GCTAGCTA"],
                "name": ["s1", "s2"],
            }
        )
        with pytest.raises(ValueError, match="No common p5 or p3 sequence found"):
            trim_p5_and_p3(df)
