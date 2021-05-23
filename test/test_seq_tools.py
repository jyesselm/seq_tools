"""
Tests for `seq_tools` module.
"""
import os
import pytest
import numpy as np
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
        "ds"       : None,
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


def test_csv_mw():
    p = get_default_args()
    base_dir = get_base_dir()
    csv_path = base_dir + "/resources/test.csv"
    p["input"] = csv_path
    p["calc"] = "mw"
    df = seq_tools.process_args(p)
    assert "molecular weight" in df.columns
    row = df.loc[1]
    assert row["molecular weight"] == 24537


def test_mult_csv():
    p = get_default_args()
    base_dir = get_base_dir()
    csv_path = base_dir + "/resources/test.csv"
    p["input"] = csv_path
    p["calc"] = "mw,len"
    df = seq_tools.process_args(p)
    assert "molecular weight" in df.columns
    assert "len" in df.columns


def test_csv_convert_to_dna():
    p = get_default_args()
    base_dir = get_base_dir()
    csv_path = base_dir + "/resources/test.csv"
    p["input"] = csv_path
    p["add_t7"] = True
    p["to_dna"] = True
    df = seq_tools.process_args(p)
    row = df.loc[1]
    assert (
            row["sequence"]
            == "TTCTAATACGACTCACTATAGATATGGATGATTAGGACATGCATTGCTGAGGGGAAACTTTTTGCAATGCAACAG"
               "CCAAATCGTCCTAAGTC"
    )


def test_trim():
    p = get_default_args()
    base_dir = get_base_dir()
    csv_path = base_dir + "/resources/test.csv"
    p["input"] = csv_path
    p["calc"] = "l"
    df = seq_tools.process_args(p)
    total = np.sum(df["len"])
    assert total == 430
    p["trim5"] = 5
    df = seq_tools.process_args(p)
    total = np.sum(df["len"])
    assert total == 400
