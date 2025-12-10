"""Comprehensive tests for CLI commands to improve coverage."""

import os
import tempfile

import pandas as pd
import pytest
from click.testing import CliRunner

from seq_tools import cli


@pytest.fixture
def temp_csv_file():
    """Create a temporary CSV file with test sequences."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("name,sequence\n")
        f.write("seq1,ATCGATCG\n")
        f.write("seq2,GCTAGCTA\n")
        f.write("seq3,TTAATTAA\n")
        temp_path = f.name
    yield temp_path
    if os.path.exists(temp_path):
        os.unlink(temp_path)


@pytest.fixture
def temp_output_file():
    """Create a temporary output file path."""
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        temp_path = f.name
    yield temp_path
    if os.path.exists(temp_path):
        os.unlink(temp_path)


class TestCLIGenerators:
    """Test CLI generator commands."""

    def test_mutate_command(self, temp_output_file):
        """Test mutate command with all options."""
        runner = CliRunner()
        result = runner.invoke(
            cli.mutate,
            [
                "ATCGATCGATCGATCG",
                "-m",
                "3",
                "-n",
                "5",
                "-nt",
                "DNA",
                "-o",
                temp_output_file,
            ],
        )
        assert result.exit_code == 0
        df = pd.read_csv(temp_output_file)
        assert len(df) == 5
        assert all(len(seq) == 16 for seq in df["sequence"])

    def test_mutate_with_primers(self, temp_output_file):
        """Test mutate command with 5' and 3' primers."""
        runner = CliRunner()
        result = runner.invoke(
            cli.mutate,
            [
                "AAAATCGATCGTTTT",
                "-m",
                "2",
                "-n",
                "3",
                "-p5",
                "AAAA",
                "-p3",
                "TTTT",
                "-o",
                temp_output_file,
            ],
        )
        assert result.exit_code == 0
        df = pd.read_csv(temp_output_file)
        assert all(seq.startswith("AAAA") for seq in df["sequence"])
        assert all(seq.endswith("TTTT") for seq in df["sequence"])

    def test_random_command(self, temp_output_file):
        """Test random sequence generation."""
        runner = CliRunner()
        result = runner.invoke(
            cli.random, ["-l", "20", "-n", "10", "-nt", "RNA", "-o", temp_output_file]
        )
        assert result.exit_code == 0
        df = pd.read_csv(temp_output_file)
        assert len(df) == 10
        assert all(len(seq) == 20 for seq in df["sequence"])
        # Check for RNA nucleotides
        assert all("U" in seq for seq in df["sequence"])

    def test_random_with_primers(self, temp_output_file):
        """Test random with constant sequences."""
        runner = CliRunner()
        result = runner.invoke(
            cli.random,
            ["-l", "30", "-n", "5", "-p5", "GGG", "-p3", "CCC", "-o", temp_output_file],
        )
        assert result.exit_code == 0
        df = pd.read_csv(temp_output_file)
        assert all(seq.startswith("GGG") for seq in df["sequence"])
        assert all(seq.endswith("CCC") for seq in df["sequence"])


class TestCLIIO:
    """Test CLI I/O commands."""

    def test_to_fasta(self, temp_csv_file):
        """Test FASTA file generation."""
        runner = CliRunner()
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as f:
            output = f.name
        try:
            result = runner.invoke(cli.to_fasta, [temp_csv_file, "-o", output])
            assert result.exit_code == 0
            assert os.path.exists(output)
            with open(output) as f:
                content = f.read()
                assert ">seq1" in content
                assert "ATCGATCG" in content
        finally:
            if os.path.exists(output):
                os.unlink(output)

    def test_to_fasta_single_sequence(self):
        """Test FASTA generation from single sequence."""
        runner = CliRunner()
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as f:
            output = f.name
        try:
            result = runner.invoke(cli.to_fasta, ["ATCGATCG", "-o", output])
            assert result.exit_code == 0
            assert os.path.exists(output)
        finally:
            if os.path.exists(output):
                os.unlink(output)


class TestCLITransforms:
    """Test CLI transform commands."""

    def test_trim_command(self, temp_csv_file, temp_output_file):
        """Test trim command."""
        runner = CliRunner()
        result = runner.invoke(
            cli.trim, [temp_csv_file, "-p5", "2", "-p3", "2", "-o", temp_output_file]
        )
        assert result.exit_code == 0
        df = pd.read_csv(temp_output_file)
        # Original sequences were 8 bases, trimming 2 from each end = 4 bases
        assert all(len(seq) == 4 for seq in df["sequence"])

    def test_fold_command(self, temp_output_file):
        """Test RNA folding command."""
        runner = CliRunner()
        result = runner.invoke(cli.fold, ["GGGGAAAACCCC", "-o", temp_output_file])
        assert result.exit_code == 0
        df = pd.read_csv(temp_output_file)
        assert "structure" in df.columns
        assert "mfe" in df.columns

    def test_transcribe_command(self, temp_output_file):
        """Test transcription command."""
        runner = CliRunner()
        # T7 promoter + some sequence with T's
        seq = "TTCTAATACGACTCACTATAGGGGTTTACCCC"
        result = runner.invoke(cli.transcribe, [seq, "-o", temp_output_file])
        assert result.exit_code == 0
        df = pd.read_csv(temp_output_file)
        # Should have removed T7 promoter (20 bases) and converted to RNA
        # Original has TTT which becomes UUU
        transcribed = df.iloc[0]["sequence"]
        assert "U" in transcribed or len(transcribed) < len(seq)
        assert "structure" in df.columns


class TestCLIAnalysis:
    """Test CLI analysis commands."""

    def test_rc_command(self, temp_output_file):
        """Test reverse complement command."""
        runner = CliRunner()
        result = runner.invoke(cli.rc, ["ATCG", "-nt", "DNA", "-o", temp_output_file])
        assert result.exit_code == 0
        df = pd.read_csv(temp_output_file)
        assert df.iloc[0]["rev_comp"] == "CGAT"

    def test_rc_rna(self, temp_output_file):
        """Test reverse complement for RNA."""
        runner = CliRunner()
        result = runner.invoke(cli.rc, ["AUCG", "-nt", "RNA", "-o", temp_output_file])
        assert result.exit_code == 0
        df = pd.read_csv(temp_output_file)
        assert df.iloc[0]["rev_comp"] == "CGAU"


class TestCLICommon:
    """Test CLI common utilities."""

    def test_get_ntype_enforcement(self, temp_csv_file, temp_output_file):
        """Test nucleotide type enforcement."""
        runner = CliRunner()
        # Test with DNA enforcement
        result = runner.invoke(cli.to_dna, [temp_csv_file, "-o", temp_output_file])
        assert result.exit_code == 0

    def test_handle_output_multiple_sequences(self, temp_csv_file):
        """Test output handling with multiple sequences."""
        runner = CliRunner()
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            output = f.name
        try:
            result = runner.invoke(cli.to_dna, [temp_csv_file, "-o", output])
            assert result.exit_code == 0
            df = pd.read_csv(output)
            assert len(df) > 1
        finally:
            if os.path.exists(output):
                os.unlink(output)


class TestCLIPrimers:
    """Test CLI primer-related commands."""

    def test_list_common_seqs_command(self):
        """Test list-common-seqs command."""
        runner = CliRunner()
        result = runner.invoke(cli.list_common_seqs)
        assert result.exit_code == 0
        # Should output tables of p5 and p3 sequences
        assert "5' sequences" in result.output or "p5" in result.output.lower()


class TestCLITransformsExtended:
    """Extended tests for transform commands."""

    def test_to_dna_template_command(self, temp_output_file):
        """Test DNA template conversion."""
        runner = CliRunner()
        result = runner.invoke(cli.to_dna_template, ["AUCG", "-o", temp_output_file])
        assert result.exit_code == 0
        df = pd.read_csv(temp_output_file)
        # Should have T7 promoter and be DNA
        assert len(df.iloc[0]["sequence"]) > 4


class TestLoggerCoverage:
    """Tests to improve logger coverage."""

    def test_logger_with_file(self):
        """Test logger with file output."""
        import tempfile

        from seq_tools.logger import setup_applevel_logger

        with tempfile.NamedTemporaryFile(suffix=".log", delete=False) as f:
            log_file = f.name
        try:
            logger = setup_applevel_logger(file_name=log_file)
            logger.info("Test message")
            assert os.path.exists(log_file)
        finally:
            if os.path.exists(log_file):
                os.unlink(log_file)

    def test_logger_debug_mode(self):
        """Test logger in debug mode."""
        import logging

        from seq_tools.logger import setup_applevel_logger

        logger = setup_applevel_logger(is_debug=True)
        assert logger.level == logging.DEBUG


class TestUtilsCoverage:
    """Tests to improve utils coverage."""

    def test_get_resources_path_fallback(self):
        """Test fallback method for resources path."""
        from seq_tools.utils import get_resources_path

        path = get_resources_path()
        assert path.exists()
        assert (path / "p5_sequences.csv").exists()


class TestDataframeIOCoverage:
    """Tests for dataframe I/O to improve coverage."""

    def test_to_opool_function(self, temp_output_file):
        """Test opool file generation directly."""
        pytest.importorskip("openpyxl")
        import pandas as pd

        from seq_tools.dataframe.io import to_opool

        df = pd.DataFrame({"name": ["seq1", "seq2"], "sequence": ["ATCG", "GCTA"]})

        output = temp_output_file.replace(".csv", ".xlsx")
        to_opool(df, "test_pool", output)
        assert os.path.exists(output)
        os.unlink(output)
