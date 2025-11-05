# seq_tools Tutorial Notebooks

This directory contains Jupyter notebooks demonstrating how to use the `seq_tools` package for working with nucleic acid sequences.

## Notebooks Overview

### 01_introduction.ipynb
- Package overview and installation
- Quick start examples
- Overview of available functions

### 02_sequence_operations.ipynb
- Working with individual sequences
- DNA/RNA conversions
- Reverse complement calculations
- Molecular weight calculations
- Extinction coefficient calculations
- Maximum stretch analysis

### 03_structure_analysis.ipynb
- RNA folding with ViennaRNA
- Working with SequenceStructure objects
- Structure-based pattern matching
- Structure-aware extinction coefficients
- Connectivity lists and pairing analysis

### 04_dataframe_operations.ipynb
- Batch processing with pandas DataFrames
- Adding/trimming sequences
- Converting between DNA and RNA
- Transcription workflows
- Batch folding
- Calculating properties for multiple sequences
- Pattern matching
- File I/O (FASTA export)

### 05_advanced_features.ipynb
- Edit distance calculations
- Parallel processing
- Generating mutated sequence libraries
- Generating random sequences
- Structure-based pattern matching in DataFrames
- Complete analysis workflows
- Data validation

## Getting Started

1. Install the package:
   ```bash
   pip install rna_seq_tools
   ```

2. Install Jupyter (if not already installed):
   ```bash
   pip install jupyter
   ```

3. Start Jupyter:
   ```bash
   jupyter notebook
   ```

4. Open the notebooks in order (01 through 05) to follow the tutorial series.

## Requirements

- Python 3.7+
- seq_tools package
- pandas
- vienna (for RNA folding)
- Jupyter notebook

All dependencies are installed automatically when you install `rna_seq_tools`.

