"""
module to test the cli.py module
"""
import os
import pandas as pd
from click.testing import CliRunner
from seq_tools import cli

resource_path = os.path.join(os.path.dirname(__file__), "resources")


def test_add():
    """
    Test the add function
    """
    runner = CliRunner()
    result = runner.invoke(cli.add, ["GGGGTTTTCCCC", "-p5", "AAAA", "-p3", "CCCC"])
    assert result.exit_code == 0
    lines = result.output.splitlines()
    assert lines[-2] == "sequence    AAAAGGGGTTTTCCCCCCCC"


def test_edit_distance():
    """
    Test the edit distance function
    """
    runner = CliRunner()
    result = runner.invoke(cli.edit_distance, f"{resource_path}/test.csv")
    assert result.exit_code == 0
    lines = result.output.splitlines()
    assert (
        lines[0] == "SEQ_TOOLS.edit_distance - INFO - edit distance: 17.666666666666668"
    )


def test_extinction_coeff():
    """
    Test the extinction_coeff function
    """
    runner = CliRunner()
    result = runner.invoke(cli.ec, ["GGGGTTTTCCCC"])
    assert result.exit_code == 0
    lines = result.output.splitlines()
    assert lines[-2] == "extinction_coeff          103700"
    result = runner.invoke(cli.ec, [f"{resource_path}/test.csv", "-nt", "RNA", "-ds"])
    assert result.exit_code == 0


def test_molecular_weight():
    """
    Test the molecular_weight function
    """
    runner = CliRunner()
    result = runner.invoke(cli.mw, ["GGGGTTTTCCCC"])
    assert result.exit_code == 0
    lines = result.output.splitlines()
    assert lines[-2] == "mw                3906.4"


def test_to_dna_single():
    """
    Test the to_dna function
    """
    runner = CliRunner()
    result = runner.invoke(cli.to_dna, ["GGGGUUUUCCCC"])
    assert result.exit_code == 0


def test_to_dna_file():
    """
    Test the to_dna function
    """
    runner = CliRunner()
    result = runner.invoke(cli.to_dna, [os.path.join(resource_path, "test.csv")])
    assert result.exit_code == 0
    assert os.path.isfile("output.csv")
    df = pd.read_csv("output.csv")
    df_expected = pd.read_csv(os.path.join(resource_path, "test_to_dna.csv"))
    os.remove("output.csv")
    assert df.equals(df_expected)


def test_to_rna_single():
    """
    Test the to_rna function
    """
    runner = CliRunner()
    result = runner.invoke(cli.to_rna, ["GGGGTTTTCCCC"])
    assert result.exit_code == 0
