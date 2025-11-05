"""
test structure module for seq_tools
"""

import pytest
from seq_tools.structure import SequenceStructure, Match, find, find_seq_struct


class TestSequenceStructure:
    """Tests for SequenceStructure class basic operations."""

    def test_init(self):
        """Test initialization of SequenceStructure object."""
        ss = SequenceStructure("ATCG", "....")
        assert ss.sequence == "ATCG"
        assert ss.structure == "...."
        # test that the sequence and structure are the same length
        with pytest.raises(ValueError):
            SequenceStructure("ATCG", "....(")

    def test_split_strands(self):
        """Test that split_strands returns a list of SequenceStructure objects."""
        ss = SequenceStructure("ATCG&ATCG", "....&....")
        assert len(ss.split_strands()) == 2
        assert isinstance(ss.split_strands()[0], SequenceStructure)

    def test_join(self):
        """Test that join returns a SequenceStructure object."""
        ss1 = SequenceStructure("ATCG", "....")
        ss2 = SequenceStructure("ATCG", "....")
        ss = ss1.join(ss2)
        assert isinstance(ss, SequenceStructure)
        assert ss.sequence == "ATCG&ATCG"
        assert ss.structure == "....&...."

    def test_replace(self):
        """Test that replace returns a SequenceStructure object."""
        ss1 = SequenceStructure("ATCG", "....")
        ss2 = SequenceStructure("AA", "()")
        ss = ss1.replace(ss2, 2)
        assert ss.sequence == "ATAA"
        assert ss.structure == "..()"


class TestFind:
    """Tests for the find function."""

    def test_find_basic(self):
        """Test that find returns the correct index."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("GAAAC", "(...)")
        r = find(struct, sub)
        assert r == [([2, 7],)]

    def test_find_multi_strand(self):
        """Test that find returns the correct index for multi-strand searches."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("GGG&CCC", "(((&)))")
        r = find(struct, sub)
        assert r == [([0, 3], [6, 9])]
        assert len(r) == 1
        # should have multiple solutions
        sub = SequenceStructure("GG&CC", "((&))")
        r = find(struct, sub)
        assert len(r) == 4

    def test_find_real_solution(self):
        """Test that find returns the correct index for a real-world example."""
        seq = (
            "GGGCUUCGGCCCACUGUCUAACGAGGAAACUUUGUUCAGAUGGAUAUUUCGUCAAUCUCGAGUAGGGA"
            "UUGAUAGAAAUAAGCGUGACGCGUCAGAGAAACUCUGACCAUCAUGUACAAAGAAACAACAACAACAAC"
        )
        ss = (
            "((((....)))).(((((((((((((...))))))).)))))).(((((((((((((((.....))))))"
            "))).)))))).((((((...(((((((...)))))))..))))))......................"
        )
        struct = SequenceStructure(seq, ss)
        sub = SequenceStructure("UCUAAC&GUUCAGA", "((((((&))).)))")
        r = find(struct, sub)
        assert r == [([16, 22], [33, 40])]
        bns = r[0]
        bmin, bmax = bns[0]
        assert struct[bmin:bmax].sequence == "UCUAAC"
        assert struct[bmin:bmax].structure == "(((((("

    def test_find_failure(self):
        """Test edge case that might fail."""
        seq = (
            "GGAACAGCACUUCGGUGCAAACAUUGAGAGCGAGUAGCUUUCAAUGAAAGCUUGUGCCCGUUGUUAUGGUUU"
            "GGGACCGAGGUUUUGAACUACUCUGAACACGGGAAACUGUACCCAGGGCGGAACCGUUUGACGUUUCGGGCCUA"
            "AGUCGGCGGGUACAGGUACAAAGAAACAACAACAACAAC"
        )
        ss = (
            "......((((....))))...((((((((((.....))))))))))...((((((((((((((....((((("
            "(((((((((((((((.....(((((...((((....))))...))))))))))))..)))..))))))))))"
            ".....))))))))))))))......................"
        )
        struct = SequenceStructure(seq, ss)
        sub = SequenceStructure("GGGAAACU", "((....))")
        r = find(struct, sub)
        print(r)


class TestFindSeqStruct:
    """Tests for the find_seq_struct function."""

    def test_alias(self):
        """Test that find_seq_struct is an alias for find and works identically."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("GAAAC", "(...)")
        r_find = find(struct, sub)
        r_find_seq_struct = find_seq_struct(struct, sub)
        assert r_find == r_find_seq_struct

    def test_no_match(self):
        """Test that find_seq_struct returns empty list when no match is found."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("TTTTT", ".....")
        r = find_seq_struct(struct, sub)
        assert r == []

    def test_multiple_matches_single_strand(self):
        """Test finding multiple occurrences on a single strand."""
        struct = SequenceStructure("GGGAAACCCGGGAAACCC", "(((...)))(((...)))")
        sub = SequenceStructure("GAAAC", "(...)")
        r = find_seq_struct(struct, sub)
        assert len(r) == 2
        # Order may vary, so check that all expected matches are present
        expected_matches = [[2, 7], [11, 16]]
        actual_matches = [match[0] for match in r]
        for expected in expected_matches:
            assert expected in actual_matches

    def test_overlapping_matches(self):
        """Test that overlapping matches are all found."""
        struct = SequenceStructure("AAAAAA", "......")
        sub = SequenceStructure("AAA", "...")
        r = find_seq_struct(struct, sub)
        # Should find 4 overlapping matches: [0,3], [1,4], [2,5], [3,6]
        assert len(r) == 4
        # Order may vary, so check that all expected matches are present
        expected_matches = [[0, 3], [1, 4], [2, 5], [3, 6]]
        actual_matches = [match[0] for match in r]
        for expected in expected_matches:
            assert expected in actual_matches

    def test_single_char_match(self):
        """Test matching a single character."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("A", ".")
        r = find_seq_struct(struct, sub)
        # Should find all three A's at positions 3, 4, 5
        assert len(r) == 3
        # Order may vary, so check that all expected matches are present
        expected_matches = [[3, 4], [4, 5], [5, 6]]
        actual_matches = [match[0] for match in r]
        for expected in expected_matches:
            assert expected in actual_matches


class TestFindSeqStructWildcards:
    """Tests for find_seq_struct with wildcard N support."""

    def test_wildcard_N(self):
        """Test that N wildcard matches any nucleotide."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        # N should match G
        sub = SequenceStructure("NGG", "(((")
        r = find_seq_struct(struct, sub)
        assert len(r) == 1
        assert r == [([0, 3],)]

        # N should match C
        sub = SequenceStructure("CCN", ")))")
        r = find_seq_struct(struct, sub)
        assert len(r) == 1
        assert r == [([6, 9],)]

    def test_multiple_N(self):
        """Test multiple N wildcards in sequence."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("NNN", "(((")
        r = find_seq_struct(struct, sub)
        assert len(r) == 1
        assert r == [([0, 3],)]


class TestFindSeqStructParameters:
    """Tests for find_seq_struct with start/end parameters."""

    def test_start_parameter(self):
        """Test that start parameter limits search range."""
        struct = SequenceStructure("GGGAAACCCGGGAAACCC", "(((...)))(((...)))")
        sub = SequenceStructure("GAAAC", "(...)")
        # Should only find second match when starting at position 10
        r = find_seq_struct(struct, sub, start=10)
        assert len(r) == 1
        assert r == [([11, 16],)]

        # Should find both when starting at 0
        r = find_seq_struct(struct, sub, start=0)
        assert len(r) == 2

    def test_end_parameter(self):
        """Test that end parameter limits search range."""
        struct = SequenceStructure("GGGAAACCCGGGAAACCC", "(((...)))(((...)))")
        sub = SequenceStructure("GAAAC", "(...)")
        # Should only find first match when ending at position 10
        r = find_seq_struct(struct, sub, end=10)
        assert len(r) == 1
        assert r == [([2, 7],)]

    def test_start_end_parameters(self):
        """Test using both start and end parameters together."""
        struct = SequenceStructure(
            "GGGAAACCCGGGAAACCCGGGAAACCC", "(((...)))(((...)))(((...)))"
        )
        sub = SequenceStructure("GAAAC", "(...)")
        # Should only find middle match
        r = find_seq_struct(struct, sub, start=5, end=18)
        assert len(r) == 1
        assert r == [([11, 16],)]

    def test_invalid_start(self):
        """Test that invalid start parameter raises ValueError."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("AAA", "...")
        with pytest.raises(ValueError, match="Invalid search range"):
            find_seq_struct(struct, sub, start=-1)
        with pytest.raises(ValueError, match="Invalid search range"):
            find_seq_struct(struct, sub, start=10)

    def test_invalid_end(self):
        """Test that invalid end parameter raises ValueError."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("AAA", "...")
        with pytest.raises(ValueError, match="Invalid search range"):
            find_seq_struct(struct, sub, end=10)
        with pytest.raises(ValueError, match="Invalid search range"):
            find_seq_struct(struct, sub, start=5, end=3)


class TestFindSeqStructBoundaries:
    """Tests for find_seq_struct boundary conditions."""

    def test_boundary_start(self):
        """Test match at the very beginning of sequence."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("GGG", "(((")
        r = find_seq_struct(struct, sub)
        assert r == [([0, 3],)]

    def test_boundary_end(self):
        """Test match at the very end of sequence."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("CCC", ")))")
        r = find_seq_struct(struct, sub)
        assert r == [([6, 9],)]

    def test_exact_match(self):
        """Test when sub-structure exactly matches full structure."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("GGGAAACCC", "(((...)))")
        r = find_seq_struct(struct, sub)
        assert r == [([0, 9],)]


class TestFindSeqStructPatterns:
    """Tests for find_seq_struct with different structure patterns."""

    def test_all_dots(self):
        """Test matching structure with all unpaired bases."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("AAA", "...")
        r = find_seq_struct(struct, sub)
        assert r == [([3, 6],)]

    def test_all_paired(self):
        """Test matching structure with all paired bases."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("GGG", "(((")
        r = find_seq_struct(struct, sub)
        assert r == [([0, 3],)]


class TestFindSeqStructMultiStrand:
    """Tests for find_seq_struct with multi-strand searches."""

    def test_multi_strand_no_match(self):
        """Test multi-strand search when no match is found."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("TTT&AAA", "(((&)))")
        r = find_seq_struct(struct, sub)
        assert r == []

    def test_multi_strand_wildcard(self):
        """Test multi-strand search with N wildcards."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("NNN&CCC", "(((&)))")
        r = find_seq_struct(struct, sub)
        assert len(r) == 1
        assert r == [([0, 3], [6, 9])]

    def test_complex_multi_strand(self):
        """Test complex multi-strand search with multiple possible combinations."""
        struct = SequenceStructure("GGGAAACCCGGGAAACCC", "(((...)))(((...)))")
        sub = SequenceStructure("GG&CC", "((&))")
        r = find_seq_struct(struct, sub)
        # Should find multiple combinations: first GG with any CC, etc.
        assert len(r) > 1
        # All matches should have format ([start, end], [start, end])
        for match in r:
            assert len(match) == 2
            assert len(match[0]) == 2
            assert len(match[1]) == 2


class TestFindSeqStructMismatches:
    """Tests for find_seq_struct with sequence/structure mismatches."""

    def test_sequence_match_structure_mismatch(self):
        """Test that sequence matches but structure mismatches return no results."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        # Sequence matches but structure doesn't (should be paired, not unpaired)
        sub = SequenceStructure("GGG", "...")
        r = find_seq_struct(struct, sub)
        assert r == []

    def test_structure_match_sequence_mismatch(self):
        """Test that structure matches but sequence mismatches return no results."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        # Structure matches but sequence doesn't
        sub = SequenceStructure("TTT", "(((")
        r = find_seq_struct(struct, sub)
        assert r == []


class TestFindSeqStructEdgeCases:
    """Tests for find_seq_struct edge cases."""

    def test_empty_sub(self):
        """Test behavior with empty sub-structure (edge case)."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("", "")
        r = find_seq_struct(struct, sub)
        # Empty string should match at all positions, but this might be implementation dependent
        # This test documents current behavior
        assert isinstance(r, list)


class TestFindSeqStructMatchFormat:
    """Tests for find_seq_struct with the new Match format."""

    def test_match_format_single_strand(self):
        """Test Match format for single-strand searches."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("GAAAC", "(...)")
        matches = find_seq_struct(struct, sub, format="match")
        assert len(matches) == 1
        assert isinstance(matches[0], Match)
        assert matches[0].start == 2
        assert matches[0].end == 7
        assert matches[0].strands == ((2, 7),)

    def test_match_format_multiple_matches(self):
        """Test Match format with multiple matches."""
        struct = SequenceStructure("GGGAAACCCGGGAAACCC", "(((...)))(((...)))")
        sub = SequenceStructure("GAAAC", "(...)")
        matches = find_seq_struct(struct, sub, format="match")
        assert len(matches) == 2
        assert all(isinstance(m, Match) for m in matches)
        # Order may vary
        starts = {m.start for m in matches}
        assert starts == {2, 11}

    def test_match_format_multi_strand(self):
        """Test Match format for multi-strand searches."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("GGG&CCC", "(((&)))")
        matches = find_seq_struct(struct, sub, format="match")
        assert len(matches) == 1
        assert isinstance(matches[0], Match)
        # Multi-strand matches don't have .start/.end attributes
        with pytest.raises(AttributeError):
            _ = matches[0].start
        assert matches[0].strands == ((0, 3), (6, 9))

    def test_match_format_no_matches(self):
        """Test Match format when no matches are found."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("TTTTT", ".....")
        matches = find_seq_struct(struct, sub, format="match")
        assert matches == []

    def test_match_format_invalid(self):
        """Test that invalid format raises ValueError."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("GAAAC", "(...)")
        with pytest.raises(ValueError, match="Invalid format"):
            find_seq_struct(struct, sub, format="invalid")

    def test_match_repr_single_strand(self):
        """Test Match string representation for single-strand."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("GAAAC", "(...)")
        matches = find_seq_struct(struct, sub, format="match")
        assert "Match(start=2, end=7)" in str(matches[0])

    def test_match_repr_multi_strand(self):
        """Test Match string representation for multi-strand."""
        struct = SequenceStructure("GGGAAACCC", "(((...)))")
        sub = SequenceStructure("GGG&CCC", "(((&)))")
        matches = find_seq_struct(struct, sub, format="match")
        assert "strands" in str(matches[0])
