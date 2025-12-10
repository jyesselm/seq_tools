"""Final push to reach 90% coverage - targeting remaining gaps."""

import os
import tempfile

import pandas as pd
import pytest

from seq_tools import dataframe, utils


class TestUtilsFinalPush:
    """Complete utils coverage."""

    def test_get_resources_path_without_files(self):
        """Test fallback when importlib.resources.files is not available."""
        import seq_tools.utils as utils_module

        # Save original value
        original = utils_module._HAS_FILES

        try:
            # Test fallback path
            utils_module._HAS_FILES = False
            path = utils.get_resources_path()
            assert path.exists()
            assert (path / "p5_sequences.csv").exists()
        finally:
            # Restore
            utils_module._HAS_FILES = original

    def test_dataframe_to_sequences(self):
        """Test dataframe to sequences conversion."""
        df = pd.DataFrame(
            {
                "name": ["seq1", "seq2"],
                "sequence": ["ATCG", "GCTA"],
                "extra": ["a", "b"],
            }
        )
        result = utils.dataframe_to_sequences(df)
        assert len(result) == 2
        assert result[0]["name"] == "seq1"
        assert result[0]["sequence"] == "ATCG"
        assert result[0]["extra"] == "a"


class TestDataframeIOFinalPush:
    """Complete dataframe I/O coverage."""

    def test_to_fasta_multiple_sequences(self):
        """Test FASTA with multiple sequences."""
        from seq_tools.dataframe.io import to_fasta

        df = pd.DataFrame(
            {"name": ["seq1", "seq2", "seq3"], "sequence": ["ATCG", "GCTA", "TTAA"]}
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            output = f.name

        try:
            to_fasta(df, output)
            with open(output) as f:
                lines = f.readlines()
                assert ">seq1\n" in lines
                assert ">seq2\n" in lines
                assert ">seq3\n" in lines
        finally:
            os.unlink(output)


class TestDataframePrimersFinalPush:
    """Complete dataframe primers coverage."""

    def test_trim_p5_and_p3_with_known_sequences(self):
        """Test trim with actual p5/p3 sequences."""
        # Create sequence with known p5/p3 from resources
        resources_path = utils.get_resources_path()
        df_p5 = pd.read_csv(resources_path / "p5_sequences.csv")
        df_p3 = pd.read_csv(resources_path / "p3_sequences.csv")

        # Use first sequences from resources
        if len(df_p5) > 0 and len(df_p3) > 0:
            p5_seq = df_p5.iloc[0]["sequence"]
            p3_seq = df_p3.iloc[0]["sequence"]

            # Create test dataframe
            df = pd.DataFrame(
                {
                    "sequence": [
                        p5_seq + "ATCGATCG" + p3_seq,
                        p5_seq + "GCTAGCTA" + p3_seq,
                    ]
                }
            )

            result = dataframe.trim_p5_and_p3(df)
            # Should have trimmed off p5 and p3
            assert len(result.iloc[0]["sequence"]) < len(df.iloc[0]["sequence"])

    def test_remove_common_p5_p3_with_structure_column(self):
        """Test structure-aware removal."""
        resources_path = utils.get_resources_path()
        df_p5 = pd.read_csv(resources_path / "p5_sequences.csv")
        df_p3 = pd.read_csv(resources_path / "p3_sequences.csv")

        # Find entries with structure info
        p5_with_struct = df_p5[df_p5["structure"].notna()]
        p3_with_struct = df_p3[df_p3["structure"].notna()]

        if len(p5_with_struct) > 0 and len(p3_with_struct) > 0:
            p5_seq = p5_with_struct.iloc[0]["sequence"]
            p5_struct = p5_with_struct.iloc[0]["structure"]
            p3_seq = p3_with_struct.iloc[0]["sequence"]
            p3_struct = p3_with_struct.iloc[0]["structure"]

            # Create test dataframe with structure
            df = pd.DataFrame(
                {
                    "sequence": [p5_seq + "ATCG" + p3_seq, p5_seq + "GCTA" + p3_seq],
                    "structure": [
                        p5_struct + "...." + p3_struct,
                        p5_struct + "...." + p3_struct,
                    ],
                }
            )

            result = dataframe.remove_common_p5_p3_by_structure(df)
            # Should have removed p5 and p3
            assert len(result.iloc[0]["sequence"]) < len(df.iloc[0]["sequence"])

    def test_remove_common_seqs_with_structure(self):
        """Test remove_common_seqs with structure info."""
        resources_path = utils.get_resources_path()
        df_p5 = pd.read_csv(resources_path / "p5_sequences.csv")
        df_p3 = pd.read_csv(resources_path / "p3_sequences.csv")

        if len(df_p5) > 0 and len(df_p3) > 0:
            p5_seq = df_p5.iloc[0]["sequence"]
            p3_seq = df_p3.iloc[0]["sequence"]

            # Create test dataframe
            df = pd.DataFrame(
                {"sequence": [p5_seq + "ATCG" + p3_seq, p5_seq + "GCTA" + p3_seq]}
            )

            result, info = dataframe.remove_common_seqs(df)
            # Should have removed p5 and p3
            assert len(result.iloc[0]["sequence"]) < len(df.iloc[0]["sequence"])
            assert "removed" in info


class TestCLIIOFinalPush:
    """Complete CLI I/O coverage."""

    def test_to_opool_with_file(self):
        """Test to_opool with CSV file input."""
        pytest.importorskip("openpyxl")
        from click.testing import CliRunner

        from seq_tools import cli

        runner = CliRunner()

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("name,sequence\n")
            f.write("seq1,ATCG\n")
            f.write("seq2,GCTA\n")
            input_file = f.name

        with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as f:
            output = f.name

        try:
            result = runner.invoke(cli.to_opool, [input_file, "-o", output])
            assert result.exit_code == 0
        finally:
            os.unlink(input_file)
            if os.path.exists(output):
                os.unlink(output)

    def test_to_fasta_with_file(self):
        """Test to_fasta with CSV file."""
        from click.testing import CliRunner

        from seq_tools import cli

        runner = CliRunner()

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("name,sequence\n")
            f.write("seq1,ATCGATCGATCG\n")
            f.write("seq2,GCTAGCTAGCTA\n")
            input_file = f.name

        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as f:
            output = f.name

        try:
            result = runner.invoke(cli.to_fasta, [input_file, "-o", output])
            assert result.exit_code == 0
            with open(output) as f:
                content = f.read()
                assert ">seq1" in content
                assert ">seq2" in content
        finally:
            os.unlink(input_file)
            os.unlink(output)
