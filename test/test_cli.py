"""
module to test the cli.py module
"""

from click.testing import CliRunner
from seq_tools import cli


def test_to_dna():
    """
    Test the to_dna function
    """
    runner = CliRunner()
    result = runner.invoke(cli.to_dna, ["GGGGUUUUCCCC"])
    assert result.exit_code == 0
