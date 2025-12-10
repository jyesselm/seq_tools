"""Edge case tests to push coverage to 90%."""

import os
import tempfile

import pandas as pd
import pytest
from click.testing import CliRunner

from seq_tools import cli


class TestDataframeIOEdgeCases:
    """Edge cases for dataframe I/O."""

    def test_to_fasta_with_multiline(self):
        """Test FASTA generation."""
        from seq_tools.dataframe.io import to_fasta

        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["A" * 100, "G" * 150],  # Long sequences
            }
        )

        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as f:
            output = f.name
        try:
            to_fasta(df, output)
            with open(output) as f:
                content = f.read()
                assert ">seq1" in content
                assert "A" * 100 in content
        finally:
            os.unlink(output)


class TestUtilsEdgeCases:
    """Edge cases for utils module."""

    def test_sequences_to_dataframe_with_names(self):
        """Test with explicit names."""
        from seq_tools.utils import sequences_to_dataframe

        sequences = ["ATCG", "GCTA", "TTAA"]
        names = ["a", "b", "c"]
        df = sequences_to_dataframe(sequences, names)
        assert list(df["name"]) == names
        assert list(df["sequence"]) == sequences

    def test_sequences_to_dataframe_mismatch_error(self):
        """Test error when sequences and names don't match."""
        from seq_tools.utils import sequences_to_dataframe

        with pytest.raises(ValueError, match="must match"):
            sequences_to_dataframe(["ATCG", "GCTA"], ["only_one_name"])


class TestDataframePrimersEdgeCases:
    """Edge cases for dataframe primers."""

    def test_has_t7_promoter_missing(self):
        """Test has_t7_promoter when missing."""
        from seq_tools.dataframe.primers import has_t7_promoter

        df = pd.DataFrame({"name": ["s1"], "sequence": ["ATCGATCG"]})
        assert not has_t7_promoter(df)

    def test_find_longest_common_prefix_empty(self):
        """Test with empty dataframe."""
        from seq_tools.dataframe.primers import find_longest_common_prefix

        df = pd.DataFrame({"sequence": []})
        assert find_longest_common_prefix(df) == ""

    def test_find_longest_common_suffix_empty(self):
        """Test with empty dataframe."""
        from seq_tools.dataframe.primers import find_longest_common_suffix

        df = pd.DataFrame({"sequence": []})
        assert find_longest_common_suffix(df) == ""

    def test_find_longest_common_prefix_single(self):
        """Test with single sequence."""
        from seq_tools.dataframe.primers import find_longest_common_prefix

        df = pd.DataFrame({"sequence": ["ATCGATCG"]})
        assert find_longest_common_prefix(df) == "ATCGATCG"

    def test_find_longest_common_suffix_single(self):
        """Test with single sequence."""
        from seq_tools.dataframe.primers import find_longest_common_suffix

        df = pd.DataFrame({"sequence": ["ATCGATCG"]})
        assert find_longest_common_suffix(df) == "ATCGATCG"

    def test_has_sequence_function(self):
        """Test has_sequence function."""
        from seq_tools.dataframe.primers import has_sequence

        df = pd.DataFrame({"sequence": ["ATCGATCG", "ATCGTTTT", "ATCGAAAA"]})
        assert has_sequence(df, "ATCG")
        assert not has_sequence(df, "GGGG")

    def test_trim_p5_and_p3(self):
        """Test trim_p5_and_p3 function."""
        from seq_tools.dataframe.primers import trim_p5_and_p3

        # This requires sequences with known p5/p3 sequences
        # For now, test it raises ValueError when no common sequences found
        df = pd.DataFrame({"sequence": ["ATCG", "GCTA"]})
        with pytest.raises(ValueError, match="No common"):
            trim_p5_and_p3(df)

    def test_remove_common_p5_p3(self):
        """Test remove_common_p5_p3 function."""
        from seq_tools.dataframe.primers import remove_common_p5_p3

        df = pd.DataFrame({"sequence": ["ATCG", "GCTA"]})
        with pytest.raises(ValueError, match="No common"):
            remove_common_p5_p3(df)

    def test_remove_common_p5_p3_by_structure(self):
        """Test remove_common_p5_p3_by_structure function."""
        from seq_tools.dataframe.primers import remove_common_p5_p3_by_structure

        df = pd.DataFrame({"sequence": ["ATCG", "GCTA"], "structure": ["....", "...."]})
        with pytest.raises(ValueError, match="No common"):
            remove_common_p5_p3_by_structure(df)

    def test_remove_common_seqs_warnings(self):
        """Test remove_common_seqs function."""
        from seq_tools.dataframe.primers import remove_common_seqs

        df = pd.DataFrame({"sequence": ["ATCG", "GCTA"], "structure": ["....", "...."]})
        with pytest.raises(ValueError, match="No common"):
            remove_common_seqs(df)

    def test_transcribe_without_t7(self):
        """Test transcribe without T7 promoter."""
        from seq_tools.dataframe.primers import transcribe

        df = pd.DataFrame({"sequence": ["ATCG"]})
        with pytest.raises(ValueError, match="T7 promoter"):
            transcribe(df, ignore_missing_t7=False)

    def test_transcribe_ignore_missing_t7(self):
        """Test transcribe ignoring missing T7."""
        from seq_tools.dataframe.primers import transcribe

        df = pd.DataFrame({"sequence": ["ATCG"]})
        result = transcribe(df, ignore_missing_t7=True)
        assert "U" in result.iloc[0]["sequence"]

    def test_add_common_seqs_no_structure(self):
        """Test add_common_seqs without structure."""
        from seq_tools.dataframe.primers import add_common_seqs

        df = pd.DataFrame({"name": ["s1"], "sequence": ["ATCG"]})
        with pytest.raises(ValueError, match="structure"):
            add_common_seqs(df, "GGG", "(((", "CCC", ")))")


class TestCLIGeneratorsEdgeCases:
    """Edge cases for CLI generators."""

    def test_mutate_from_file(self):
        """Test mutate command from file."""
        runner = CliRunner()
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("name,sequence\nseq1,ATCGATCGATCGATCG\n")
            input_file = f.name

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            output = f.name

        try:
            result = runner.invoke(
                cli.mutate, [input_file, "-m", "2", "-n", "3", "-o", output]
            )
            assert result.exit_code == 0
        finally:
            os.unlink(input_file)
            os.unlink(output)

    def test_mutate_error_cases(self):
        """Test mutate error handling."""
        runner = CliRunner()
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("name,sequence\nseq1,ATG\nseq2,GGG\n")
            input_file = f.name

        try:
            # Should fail - multiple sequences in template file
            result = runner.invoke(cli.mutate, [input_file, "-m", "2", "-n", "3"])
            assert result.exit_code != 0
        finally:
            os.unlink(input_file)


class TestCLIIOEdgeCases:
    """Edge cases for CLI I/O commands."""

    def test_to_opool_command(self):
        """Test to_opool command."""
        pytest.importorskip("openpyxl")
        runner = CliRunner()

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("name,sequence\nseq1,ATCG\nseq2,GCTA\n")
            input_file = f.name

        with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as f:
            output = f.name

        try:
            result = runner.invoke(
                cli.to_opool, [input_file, "-n", "test_pool", "-o", output]
            )
            assert result.exit_code == 0
            assert os.path.exists(output)
        finally:
            os.unlink(input_file)
            if os.path.exists(output):
                os.unlink(output)


class TestDataframeGeneratorsEdgeCases:
    """Edge cases for dataframe generators."""

    def test_generate_mutated_sequences_error_too_short(self):
        """Test error when template is too short."""
        from seq_tools.dataframe.generators import generate_mutated_sequences

        with pytest.raises(ValueError, match="shorter"):
            generate_mutated_sequences("ATCG", 10, 5, p5_seq="AA", p3_seq="TT")

    def test_generate_random_sequences_error_too_short(self):
        """Test error when length is too short for primers."""
        from seq_tools.dataframe.generators import generate_random_sequences

        with pytest.raises(ValueError, match="must be greater"):
            generate_random_sequences(5, 10, p5_seq="AAAA", p3_seq="TTTT")
