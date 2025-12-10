"""
Command-line interface for seq_tools.

This package provides a CLI for manipulating sequences in CSV files.
The CLI is organized into modules by functionality:

- common: Shared utilities for all commands
- transforms: Transform commands (add, fold, to_dna, to_rna, etc.)
- analysis: Analysis commands (ec, mw, rc, edit_distance)
- primers: Primer commands (has_p5, has_p3, add/remove/list common seqs)
- generators: Generator commands (mutate, random)
- io: I/O commands (to_fasta, to_opool)
"""

import click

from seq_tools.cli import analysis, generators, io, primers, transforms


@click.group()
def cli() -> None:
    """A set of scripts to manipulate sequences in CSV files."""


# Register transform commands
cli.add_command(transforms.add)
cli.add_command(transforms.fold)
cli.add_command(transforms.to_dna)
cli.add_command(transforms.to_dna_template)
cli.add_command(transforms.to_rna)
cli.add_command(transforms.trim)
cli.add_command(transforms.transcribe)

# Register analysis commands
cli.add_command(analysis.edit_distance)
cli.add_command(analysis.ec)
cli.add_command(analysis.mw)
cli.add_command(analysis.rc)

# Register primer commands
cli.add_command(primers.has_p5)
cli.add_command(primers.has_p3)
cli.add_command(primers.add_common_seqs)
cli.add_command(primers.list_common_seqs)
cli.add_command(primers.remove_common_seqs)

# Register generator commands
cli.add_command(generators.mutate)
cli.add_command(generators.random)

# Register I/O commands
cli.add_command(io.to_fasta)
cli.add_command(io.to_opool)


# Export all commands for backward compatibility
__all__ = [
    "cli",
    "add",
    "fold",
    "to_dna",
    "to_dna_template",
    "to_rna",
    "trim",
    "transcribe",
    "edit_distance",
    "ec",
    "mw",
    "rc",
    "has_p5",
    "has_p3",
    "add_common_seqs",
    "list_common_seqs",
    "remove_common_seqs",
    "mutate",
    "random",
    "to_fasta",
    "to_opool",
]

# Import commands for backward compatibility
add = transforms.add
fold = transforms.fold
to_dna = transforms.to_dna
to_dna_template = transforms.to_dna_template
to_rna = transforms.to_rna
trim = transforms.trim
transcribe = transforms.transcribe

edit_distance = analysis.edit_distance
ec = analysis.ec
mw = analysis.mw
rc = analysis.rc

has_p5 = primers.has_p5
has_p3 = primers.has_p3
add_common_seqs = primers.add_common_seqs
list_common_seqs = primers.list_common_seqs
remove_common_seqs = primers.remove_common_seqs

mutate = generators.mutate
random = generators.random

to_fasta = io.to_fasta
to_opool = io.to_opool


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    cli()
