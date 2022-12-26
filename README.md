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
Adds a sequence to the 5' and/or 3' end of a sequence. 
```shell
$ seq_tools add -p5 "AAAA" "GGGGUUUUCCCC"
SEQ_TOOLS.get_input_dataframe - INFO - reading sequence GGGGUUUUCCCC
SEQ_TOOLS.handle_output - INFO - output->
name                     seq
sequence    AAAAGGGGUUUUCCCC
Name: 0, dtype: object
```

### ec 
Calculate the extinction coefficient for each sequence. 
```shell
$ seq-tools ec "GGGGUUUUCCCC"
SEQ_TOOLS.get_input_dataframe - INFO - reading sequence GGGGUUUUCCCC
SEQ_TOOLS.handle_ntype - INFO - determining nucleic acid type: RNA
SEQ_TOOLS.handle_output - INFO - output->
name                         seq
sequence            GGGGUUUUCCCC
extinction_coeff          109500
Name: 0, dtype: object
```

### edit-distance
Calculate the edit distance of a library. On average how different each sequence 
is from the rest of the library. 
```shell
seq-tools edit-distance test/resources/test.csv
SEQ_TOOLS.edit_distance - INFO - edit distance: 17.666666666666668
```

### fold
Fold rna sequences. 
```shell
$ seq-tools fold "GGGGUUUUCCCC"
SEQ_TOOLS.get_input_dataframe - INFO - reading sequence GGGGUUUUCCCC
SEQ_TOOLS.handle_output - INFO - output->
name                   seq
sequence      GGGGUUUUCCCC
structure     ((((....))))
mfe                   -5.9
ens_defect            0.38
Name: 0, dtype: object
```

### to-dna
Convert all sequences to DNA i.e. replace T with U. 
```shell
$ seq_tools to-dna "GGGGUUUUCCCC"
SEQ_TOOLS.get_input_dataframe - INFO - reading sequence GGGGUUUUCCCC
SEQ_TOOLS.to_dna - INFO - converted sequence: GGGGTTTTCCCC
```

### other non commandline

#### structure representation

```python
from seq_tools import SequenceStructure
struct = SequenceStructure("GGGGUUUUCCCC", "((((....))))")
```
