"""
module to test the cli.py module
"""
import os
import pandas as pd
from click.testing import CliRunner
from seq_tools import cli

resource_path = os.path.join(os.path.dirname(__file__), "resources")


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
    result = runner.invoke(
        cli.to_dna, [os.path.join(resource_path, "test.csv")]
    )
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
