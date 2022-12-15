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

```shell
$ seq_tools to-dna "GGGGUUUUCCCC"
SEQ_TOOLS.get_input_dataframe - INFO - reading sequence GGGGUUUUCCCC
SEQ_TOOLS.to_dna - INFO - converted sequence: GGGGTTTTCCCC
```