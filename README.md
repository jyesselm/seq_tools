# seq_tools

[![PyPI version](https://badge.fury.io/py/rna_seq_tools.svg)](https://badge.fury.io/py/rna_seq_tools)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/jyesselm/seq_tools/actions/workflows/tests.yml/badge.svg)](https://github.com/jyesselm/seq_tools/actions)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![License](https://img.shields.io/badge/license-Non--Commercial-red.svg)](LICENSE)

A Python package for manipulating and analyzing nucleic acid sequences (DNA and RNA) in pandas DataFrames.

## Features

- **Batch operations**: Work with sequences in pandas DataFrames for efficient processing
- **Sequence manipulation**: Convert between DNA/RNA, reverse complement, add sequences
- **Structure prediction**: Fold RNA sequences using ViennaRNA
- **Analysis tools**: Calculate molecular weights, extinction coefficients, edit distances
- **CLI interface**: Command-line tools for quick sequence operations
- **Python API**: Full programmatic access to all functionality

## Installation

```bash
pip install rna_seq_tools
```

## Quick Start

### Command Line Interface

```bash
# Get help
seq-tools --help

# Convert RNA to DNA
seq-tools to-dna "AUCG"

# Fold RNA sequence
seq-tools fold "GGGGUUUUCCCC"

# Calculate molecular weight
seq-tools mw "ATCG"
```

### Python API

```python
import pandas as pd
from seq_tools import sequences_to_dataframe, fold, get_molecular_weight_df, to_rna_df

# Create a DataFrame from sequences
sequences = ["ATCG", "GCTA", "AAAA"]
df = sequences_to_dataframe(sequences)

# Convert to RNA
df = to_rna_df(df)

# Fold RNA sequences
df = fold(df)

# Calculate molecular weights
df = get_molecular_weight_df(df, "RNA", double_stranded=False)

print(df)
```

### Single Sequence Functions

For single sequence operations, import from the `sequence` module:

```python
from seq_tools.sequence import to_dna, to_rna, get_reverse_complement, get_molecular_weight

# Convert sequences
rna_seq = to_rna("ATCG")  # Returns "AUCG"
dna_seq = to_dna("AUCG")  # Returns "ATCG"

# Reverse complement
rc = get_reverse_complement("ATCG", "DNA")  # Returns "CGAT"

# Molecular weight
mw = get_molecular_weight("ATCG", "DNA")  # Returns 1307.80
```

## CLI Commands

### `add`
Add a sequence to the 5' and/or 3' end of sequences.

```bash
seq-tools add -p5 "AAAA" "GGGGUUUUCCCC"
seq-tools add -p5 "AAAA" -p3 "CCCC" input.csv
```

### `ec`
Calculate the extinction coefficient for each sequence.

```bash
seq-tools ec "GGGGUUUUCCCC"
seq-tools ec input.csv -nt RNA -ds  # RNA, double-stranded
```

### `edit-distance`
Calculate the average edit distance of a sequence library.

```bash
seq-tools edit-distance input.csv
seq-tools edit-distance input.csv --parallel --workers 4
```

### `fold`
Fold RNA sequences using ViennaRNA.

```bash
seq-tools fold "GGGGUUUUCCCC"
seq-tools fold input.csv
```

### `mw`
Calculate the molecular weight for each sequence.

```bash
seq-tools mw "ATCG"
seq-tools mw input.csv -nt DNA -ds  # DNA, double-stranded
```

### `rc`
Calculate reverse complement for each sequence.

```bash
seq-tools rc "ATCG"
seq-tools rc input.csv -nt DNA
```

### `to-dna`
Convert RNA sequences to DNA (replace U with T).

```bash
seq-tools to-dna "AUCG"
seq-tools to-dna input.csv -o output.csv
```

### `to-dna-template`
Convert RNA sequences to DNA template with T7 promoter.

```bash
seq-tools to-dna-template "AUCG"
seq-tools to-dna-template input.csv
```

### `to-rna`
Convert DNA sequences to RNA (replace T with U).

```bash
seq-tools to-rna "ATCG"
seq-tools to-rna input.csv
```

### `transcribe`
Transcribe DNA template sequences to RNA (removes T7 promoter).

```bash
seq-tools transcribe input.csv
```

### `trim`
Trim 5'/3' ends of sequences.

```bash
seq-tools trim input.csv --start 5 --end 3
```

### `to-fasta`
Generate FASTA file from CSV.

```bash
seq-tools to-fasta input.csv output.fasta
```

### `to-opool`
Generate oligo pool file (Excel) from CSV.

```bash
seq-tools to-opool input.csv "pool_name" output.xlsx
```

## DataFrame Functions

The package provides comprehensive DataFrame operations:

- **Conversion**: `to_dna_df()`, `to_rna_df()`, `to_dna_template_df()`
- **Analysis**: `get_molecular_weight_df()`, `get_extinction_coeff()`, `get_length()`
- **Structure**: `fold()` - predict RNA secondary structures
- **Manipulation**: `add()`, `trim()`, `get_reverse_complement_df()`
- **Generation**: `generate_random_sequences()`, `generate_mutated_sequences()`
- **Validation**: `has_t7_promoter()`, `has_5p_sequence()`, `has_3p_sequence()`
- **File I/O**: `to_fasta()`, `to_opool()`

See the [notebooks](notebooks/) directory for detailed examples.

## Requirements

- Python 3.7+
- pandas
- numpy
- ViennaRNA (for structure prediction)
- editdistance
- click
- tabulate

## Tutorial Notebooks

Interactive Jupyter notebooks are available in the [`notebooks/`](notebooks/) directory:

- **01_introduction.ipynb**: Package overview and quick start
- **02_sequence_operations.ipynb**: Working with individual sequences
- **03_structure_analysis.ipynb**: RNA folding and structure analysis
- **04_dataframe_operations.ipynb**: Batch processing with DataFrames
- **05_advanced_features.ipynb**: Advanced features and workflows

See the [notebooks README](notebooks/README.md) for more details.

## Development

```bash
# Clone the repository
git clone https://github.com/jyesselm/seq_tools.git
cd seq_tools

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in editable mode
pip install -e .

# Run tests
pytest test/ -v
```

## License

This project is licensed under a **Non-Commercial License**. Commercial use is prohibited. See [LICENSE](LICENSE) file for details.

For commercial licensing inquiries, please contact jyesselm@unl.edu.

## Author

**Joe Yesselman** - jyesselm@unl.edu

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
