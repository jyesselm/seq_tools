"""
simple functions for manipulating sequences and secondary structures
"""

import re
import itertools
from dataclasses import dataclass
import pandas as pd


@dataclass(frozen=True, order=True)
class SequenceStructure:
    """
    A class to hold the parameters for a structure
    """

    sequence: str
    structure: str

    def __post_init__(self):
        """
        check that the sequence and structure are the same length
        """
        if len(self.sequence) != len(self.structure):
            raise ValueError(
                f"sequence and structure are not the same length:"
                f" {self.sequence} {self.structure}"
            )

    def __add__(self, other):
        """
        add two StructParams objects
        """
        return SequenceStructure(
            self.sequence + other.sequence, self.structure + other.structure
        )

    def __len__(self):
        """
        return the length of the sequence
        """
        return len(self.sequence)

    def __getitem__(self, item):
        """
        return the attribute of the sequence
        """
        return SequenceStructure(self.sequence[item], self.structure[item])

    def to_dict(self):
        """
        return a dictionary representation of the object
        """
        return {"sequence": self.sequence, "structure": self.structure}

    def to_comma_deliminted(self):
        """
        return a csv representation of the object
        """
        return f"{self.sequence},{self.structure}"

    def split_strands(self):
        """
        split both sequence and structure over `&` and return a list of
        of StructParams objects
        :return: list of strands
        """
        seqs = self.sequence.split("&")
        structs = self.structure.split("&")
        return [SequenceStructure(s, st) for s, st in zip(seqs, structs)]

    def insert(self, pos, other):
        """
        insert a StructParams object at a given position
        """
        seq = self.sequence[:pos] + other.sequence + self.sequence[pos:]
        struct = self.structure[:pos] + other.structure + self.structure[pos:]
        return SequenceStructure(seq, struct)

    def join(self, other):
        """
        join two StructParams objects with a "&" seperating each strand
        """
        return SequenceStructure(
            self.sequence + "&" + other.sequence,
            self.structure + "&" + other.structure,
        )

    def replace(self, other, pos):
        """
        Replaces the sequence and structure of this object with the sequence and
        structure in the supplied SequenceStructure object at the specified position.

        :param other: The SequenceStructure object containing the new sequence and
        structure.
        :type other: SequenceStructure
        :param pos: The position at which to replace the sequence and structure.
        :type pos: int
        """
        if pos < 0 or pos > len(self.sequence):
            raise ValueError(f"Invalid position: {pos}")
        sequence = (
            self.sequence[:pos]
            + other.sequence
            + self.sequence[pos + len(other.sequence) :]
        )
        structure = (
            self.structure[:pos]
            + other.structure
            + self.structure[pos + len(other.structure) :]
        )
        return SequenceStructure(sequence, structure)


def find(struct: SequenceStructure, sub: SequenceStructure, start=None, end=None):
    """
    find the position of a substructure in a structure
    :param struct: the structure to search
    :param sub: the substructure to search for
    :param start: the start position to search from
    :param end: the end position to search to
    """
    if start is None:
        start = 0
    if end is None:
        end = len(struct)
    struct = struct[start:end]
    strands = sub.split_strands()
    strand_matches = []
    for strand in strands:
        pattern_seq = re.compile(
            (r"(?=(" + strand.sequence.replace("N", r"\S") + r"))")
        )
        pattern_ss = re.compile(
            (
                r"(?=("
                + strand.structure.replace("(", r"\(")
                .replace(")", r"\)")
                .replace(".", r"\.")
                + r"))"
            )
        )
        matches_seq = [
            str(m.start() + start) + "-" + str(m.end() + len(m.group(1)) + start)
            for m in pattern_seq.finditer(struct.sequence)
        ]
        matches_ss = [
            str(m.start() + start) + "-" + str(m.end() + len(m.group(1)) + start)
            for m in pattern_ss.finditer(struct.structure)
        ]
        matches = list(set(matches_seq).intersection(set(matches_ss)))
        # split each string in inter into a list of ints
        matches = [list(map(int, i.split("-"))) for i in matches]
        strand_matches.append(matches)
    all_matches = list(itertools.product(*strand_matches))
    return all_matches
