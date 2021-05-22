"""
Tests for `seq_tools` module.
"""
import os
import pytest
from seq_tools import seq_tools


def get_base_dir():
    file_path = os.path.realpath(__file__)
    spl = file_path.split("/")
    base_dir = "/".join(spl[:-1])
    return base_dir


def get_default_args():
    p = {
        "ss"       : None,
        "trim5"    : None,
        "trim3"    : None,
        "add5"     : None,
        "add3"     : None,
        "output"   : "output.csv",
        "to_rna"   : False,
        "to_dna"   : False,
        "remove_t7": False,
        "add_t7"   : False,
        "ds"       : False,
        "fold"     : False,
        "calc"     : None
    }
    return p


# tests ##############################################################################

def test_seq_len():
    p = get_default_args()
    p["input"] = "GGGAAACCC"
    p["calc"] = "len"
    df = seq_tools.process_args(p)
    row = df.loc[1]
    assert "len" in df.columns
    assert row["len"] == 9