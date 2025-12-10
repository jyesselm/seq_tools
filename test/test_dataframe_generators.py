"""
Tests for dataframe.generators module
"""

import pytest

from seq_tools.dataframe.generators import (
    generate_mutated_sequences,
    generate_random_sequences,
)


class TestGenerateMutatedSequences:
    """Tests for generate_mutated_sequences function."""

    def test_generate_mutated_sequences_basic(self):
        """Test basic mutation generation."""
        template = "ATCGATCGATCG"
        num_mutations = 2
        num_sequences = 5

        df = generate_mutated_sequences(template, num_mutations, num_sequences)

        assert len(df) == num_sequences
        assert "name" in df.columns
        assert "sequence" in df.columns

        # Check that all sequences have the same length as template
        for seq in df["sequence"]:
            assert len(seq) == len(template)

    def test_generate_mutated_sequences_names(self):
        """Test that generated names follow expected pattern."""
        template = "ATCGATCG"
        df = generate_mutated_sequences(template, 1, 3)

        expected_names = ["mutated_seq_1", "mutated_seq_2", "mutated_seq_3"]
        assert list(df["name"]) == expected_names

    def test_generate_mutated_sequences_dna(self):
        """Test mutation with DNA nucleotides."""
        template = "AAAAAAAAAA"
        num_mutations = 5
        df = generate_mutated_sequences(template, num_mutations, 10, ntype="DNA")

        # All sequences should only contain DNA nucleotides
        for seq in df["sequence"]:
            assert all(nt in "ATCG" for nt in seq)
            # Should have mutations (not identical to template)
            assert seq != template

    def test_generate_mutated_sequences_rna(self):
        """Test mutation with RNA nucleotides."""
        template = "AAAAAAAAAAUUUUUUUUUU"
        num_mutations = 3
        df = generate_mutated_sequences(template, num_mutations, 5, ntype="RNA")

        # All sequences should only contain RNA nucleotides
        for seq in df["sequence"]:
            assert all(nt in "AUCG" for nt in seq)

    def test_generate_mutated_sequences_with_p5_seq(self):
        """Test mutation with constant 5' sequence."""
        template = "GGGGGGGGGGGGGG"
        p5_seq = "AAAA"
        num_mutations = 2
        df = generate_mutated_sequences(
            template, num_mutations, 5, p5_seq=p5_seq, ntype="DNA"
        )

        # All sequences should start with p5_seq
        for seq in df["sequence"]:
            assert seq.startswith(p5_seq)
            assert len(seq) == len(template)

    def test_generate_mutated_sequences_with_p3_seq(self):
        """Test mutation with constant 3' sequence."""
        template = "GGGGGGGGGGGGGG"
        p3_seq = "TTTT"
        num_mutations = 2
        df = generate_mutated_sequences(
            template, num_mutations, 5, p3_seq=p3_seq, ntype="DNA"
        )

        # All sequences should end with p3_seq
        for seq in df["sequence"]:
            assert seq.endswith(p3_seq)
            assert len(seq) == len(template)

    def test_generate_mutated_sequences_with_both_p5_p3(self):
        """Test mutation with both constant 5' and 3' sequences."""
        template = "AAAAAAAAGGGGGGGG"
        p5_seq = "AAAA"
        p3_seq = "GGGG"
        num_mutations = 2
        df = generate_mutated_sequences(
            template, num_mutations, 5, p5_seq=p5_seq, p3_seq=p3_seq, ntype="DNA"
        )

        # All sequences should have constant ends
        for seq in df["sequence"]:
            assert seq.startswith(p5_seq)
            assert seq.endswith(p3_seq)
            assert len(seq) == len(template)

    def test_generate_mutated_sequences_variable_region_too_short(self):
        """Test error when variable region is shorter than num_mutations."""
        template = "AAAAAATTTT"
        p5_seq = "AAAA"
        p3_seq = "TTTT"
        # Variable region = 2, but num_mutations = 5
        with pytest.raises(ValueError, match="Variable region.*is shorter than"):
            generate_mutated_sequences(
                template, 5, 3, p5_seq=p5_seq, p3_seq=p3_seq, ntype="DNA"
            )

    def test_generate_mutated_sequences_template_too_short_for_p5_p3(self):
        """Test error when template is shorter than combined p5 and p3."""
        template = "ATCG"
        p5_seq = "AAAAA"
        p3_seq = "TTTTT"
        with pytest.raises(ValueError, match="Template sequence.*is shorter than"):
            generate_mutated_sequences(
                template, 1, 3, p5_seq=p5_seq, p3_seq=p3_seq, ntype="DNA"
            )

    def test_generate_mutated_sequences_all_positions_mutated(self):
        """Test mutation when num_mutations equals variable region length."""
        template = "AAAA"
        num_mutations = 4
        df = generate_mutated_sequences(template, num_mutations, 5, ntype="DNA")

        # All positions should be mutated, so no sequence should be all A's
        for seq in df["sequence"]:
            assert seq != template

    def test_generate_mutated_sequences_single_mutation(self):
        """Test generation with single mutation."""
        template = "AAAAAAAAAA"
        df = generate_mutated_sequences(template, 1, 10, ntype="DNA")

        for seq in df["sequence"]:
            # Count differences from template
            differences = sum(1 for i in range(len(seq)) if seq[i] != template[i])
            assert differences == 1

    def test_generate_mutated_sequences_zero_mutations(self):
        """Test generation with zero mutations (edge case)."""
        template = "ATCGATCG"
        df = generate_mutated_sequences(template, 0, 3, ntype="DNA")

        # All sequences should be identical to template
        for seq in df["sequence"]:
            assert seq == template


class TestGenerateRandomSequences:
    """Tests for generate_random_sequences function."""

    def test_generate_random_sequences_basic(self):
        """Test basic random sequence generation."""
        length = 20
        num_sequences = 5

        df = generate_random_sequences(length, num_sequences)

        assert len(df) == num_sequences
        assert "name" in df.columns
        assert "sequence" in df.columns

        # Check that all sequences have the correct length
        for seq in df["sequence"]:
            assert len(seq) == length

    def test_generate_random_sequences_names(self):
        """Test that generated names follow expected pattern."""
        df = generate_random_sequences(10, 3)

        expected_names = ["random_seq_1", "random_seq_2", "random_seq_3"]
        assert list(df["name"]) == expected_names

    def test_generate_random_sequences_dna(self):
        """Test random sequence generation with DNA."""
        df = generate_random_sequences(50, 10, ntype="DNA")

        # All sequences should only contain DNA nucleotides
        for seq in df["sequence"]:
            assert all(nt in "ATCG" for nt in seq)

    def test_generate_random_sequences_rna(self):
        """Test random sequence generation with RNA."""
        df = generate_random_sequences(50, 10, ntype="RNA")

        # All sequences should only contain RNA nucleotides
        for seq in df["sequence"]:
            assert all(nt in "AUCG" for nt in seq)

    def test_generate_random_sequences_with_p5_seq(self):
        """Test random generation with constant 5' sequence."""
        length = 20
        p5_seq = "AAAA"
        df = generate_random_sequences(length, 5, p5_seq=p5_seq, ntype="DNA")

        # All sequences should start with p5_seq and have correct total length
        for seq in df["sequence"]:
            assert seq.startswith(p5_seq)
            assert len(seq) == length

    def test_generate_random_sequences_with_p3_seq(self):
        """Test random generation with constant 3' sequence."""
        length = 20
        p3_seq = "TTTT"
        df = generate_random_sequences(length, 5, p3_seq=p3_seq, ntype="DNA")

        # All sequences should end with p3_seq and have correct total length
        for seq in df["sequence"]:
            assert seq.endswith(p3_seq)
            assert len(seq) == length

    def test_generate_random_sequences_with_both_p5_p3(self):
        """Test random generation with both constant 5' and 3' sequences."""
        length = 20
        p5_seq = "AAAA"
        p3_seq = "TTTT"
        df = generate_random_sequences(
            length, 5, p5_seq=p5_seq, p3_seq=p3_seq, ntype="DNA"
        )

        # All sequences should have constant ends and correct total length
        for seq in df["sequence"]:
            assert seq.startswith(p5_seq)
            assert seq.endswith(p3_seq)
            assert len(seq) == length

    def test_generate_random_sequences_length_too_short(self):
        """Test error when total length is not greater than p5 + p3 lengths."""
        length = 10
        p5_seq = "AAAAAAA"  # 7 nt
        p3_seq = "TTTT"  # 4 nt
        # Total = 11 which is > 10
        with pytest.raises(
            ValueError, match="Sequence length.*must be greater than the sum"
        ):
            generate_random_sequences(length, 3, p5_seq=p5_seq, p3_seq=p3_seq)

    def test_generate_random_sequences_length_equals_p5_p3(self):
        """Test error when length exactly equals p5 + p3."""
        length = 8
        p5_seq = "AAAA"
        p3_seq = "TTTT"
        with pytest.raises(
            ValueError, match="Sequence length.*must be greater than the sum"
        ):
            generate_random_sequences(length, 3, p5_seq=p5_seq, p3_seq=p3_seq)

    def test_generate_random_sequences_single_nucleotide_random_region(self):
        """Test generation when random region is only 1 nucleotide."""
        length = 9
        p5_seq = "AAAA"
        p3_seq = "TTTT"
        df = generate_random_sequences(
            length, 5, p5_seq=p5_seq, p3_seq=p3_seq, ntype="DNA"
        )

        for seq in df["sequence"]:
            assert len(seq) == length
            assert seq[:4] == p5_seq
            assert seq[-4:] == p3_seq
            # Middle should be 1 random nucleotide
            assert seq[4] in "ATCG"

    def test_generate_random_sequences_short_length(self):
        """Test generation of very short sequences."""
        df = generate_random_sequences(5, 10, ntype="DNA")

        for seq in df["sequence"]:
            assert len(seq) == 5
            assert all(nt in "ATCG" for nt in seq)

    def test_generate_random_sequences_long_length(self):
        """Test generation of long sequences."""
        length = 1000
        df = generate_random_sequences(length, 3, ntype="DNA")

        for seq in df["sequence"]:
            assert len(seq) == length
            assert all(nt in "ATCG" for nt in seq)

    def test_generate_random_sequences_single_sequence(self):
        """Test generation of a single sequence."""
        df = generate_random_sequences(20, 1, ntype="DNA")

        assert len(df) == 1
        assert df.iloc[0]["name"] == "random_seq_1"
        assert len(df.iloc[0]["sequence"]) == 20

    def test_generate_random_sequences_are_different(self):
        """Test that generated sequences are actually random (different)."""
        df = generate_random_sequences(50, 20, ntype="DNA")

        # Check that we have some diversity in sequences
        sequences = list(df["sequence"])
        unique_sequences = set(sequences)
        # With 20 random 50-nt sequences, we should have mostly unique sequences
        assert len(unique_sequences) > 15
