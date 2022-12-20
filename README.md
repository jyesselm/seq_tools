# seq_tools

[![PYPI package](https://badge.fury.io/py/rna_seq_tools.png)](http://badge.fury.io/py/rna_seq_tools)
[![linting: pylint](https://img.shields.io/badge/linting-pylint-yellowgreen)](https://github.com/PyCQA/pylint)
[![formatting: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

a short python tool for working with sequences in dataframes

## how to install

```shell
pip install rna_seq_tools
```

## how to use

`seq_tools` is a python package that contains a few functions for working with sequences in
dataframes. If there is a single sequence results are printed. If input is a csv then a new csv is
created with the results. Default output is "output.csv" but can be changed with the `-o` flag.

```shell
$ seq_tools --help
Usage: seq_tools [OPTIONS] COMMAND [ARGS]...

  a set scripts to manipulate sequences in csv files

Options:
  --help  Show this message and exit.

Commands:
  add              add a sequence to 5' and/or 3'
  ec               calculate the extinction coefficient for each sequence
  edit-distance    calculate the edit distance of a library
  fold             fold rna sequences
  mw               calculate the molecular weight for each sequence
  rc               calculate reverse complement for each sequence
  to-dna           convert rna sequence(s) to dna
  to-dna-template  convert rna sequence(s) to dna template, includes T7...
  to-fasta         generate fasta file from csv
  to-opool         generate oligo pool file from csv
  to-rna           convert rna sequence(s) to dna
  transcribe       convert dna sequence(s) to rna
  trim             trim 5'/3' ends of sequences

```

### add

```shell
$ seq_tools add -p5 "AAAA" "GGGGUUUUCCCC"
SEQ_TOOLS.get_input_dataframe - INFO - reading sequence GGGGUUUUCCCC
SEQ_TOOLS.handle_output - INFO - output->
name                     seq
sequence    AAAAGGGGUUUUCCCC
Name: 0, dtype: object
```

```shell
$ seq_tools to-dna "GGGGUUUUCCCC"
SEQ_TOOLS.get_input_dataframe - INFO - reading sequence GGGGUUUUCCCC
SEQ_TOOLS.to_dna - INFO - converted sequence: GGGGTTTTCCCC
```