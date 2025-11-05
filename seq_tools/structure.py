"""Simple functions for manipulating sequences and secondary structures."""

import re
import itertools
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from seq_tools.sequence import get_reverse_complement

# Constant for strand separator in multi-strand structures
STRAND_SEPARATOR = "&"


@dataclass(frozen=True, order=True)
class Match:
    """Represents a match found by find_seq_struct or find.

    A Match represents one occurrence of a substructure within a structure.
    For single-strand searches, strands contains one (start, end) tuple.
    For multi-strand searches, strands contains multiple (start, end) tuples,
    one for each strand in the match.

    Attributes:
        strands: List of (start, end) tuples representing the position of each
                 strand in the match. For single-strand matches, contains one tuple.

    Examples:
        >>> struct = SequenceStructure("GGGAAACCC", "(((...)))")
        >>> sub = SequenceStructure("GAAAC", "(...)")
        >>> matches = find_seq_struct(struct, sub, format="match")
        >>> match = matches[0]
        >>> match.start  # For single-strand matches
        2
        >>> match.end
        7
        >>> match.strands  # Can also access as list
        [(2, 7)]
    """

    strands: tuple[tuple[int, int], ...]

    def __post_init__(self) -> None:
        """Validate that strands is not empty."""
        if not self.strands:
            raise ValueError("Match must contain at least one strand")

    @property
    def start(self) -> int:
        """Start position of the match (for single-strand matches).

        Returns:
            Start position of the first strand.

        Raises:
            AttributeError: If this is a multi-strand match (use .strands instead).
        """
        if len(self.strands) > 1:
            raise AttributeError(
                "Cannot use .start on multi-strand match. Use .strands instead."
            )
        return self.strands[0][0]

    @property
    def end(self) -> int:
        """End position of the match (for single-strand matches).

        Returns:
            End position of the first strand.

        Raises:
            AttributeError: If this is a multi-strand match (use .strands instead).
        """
        if len(self.strands) > 1:
            raise AttributeError(
                "Cannot use .end on multi-strand match. Use .strands instead."
            )
        return self.strands[0][1]

    def __repr__(self) -> str:
        """Return string representation of the match."""
        if len(self.strands) == 1:
            start, end = self.strands[0]
            return f"Match(start={start}, end={end})"
        return f"Match(strands={self.strands})"


@dataclass(frozen=True, order=True)
class SequenceStructure:
    """Represents a nucleic acid sequence paired with its secondary structure.

    This class holds a sequence (DNA/RNA) and its corresponding dot-bracket
    structure notation. It supports various operations like slicing, concatenation,
    and searching for substructures.

    Attributes:
        sequence: The nucleotide sequence (DNA or RNA).
        structure: The secondary structure in dot-bracket notation (e.g., "(((...)))").

    Examples:
        >>> ss = SequenceStructure("GGGAAACCC", "(((...)))")
        >>> len(ss)
        9
        >>> ss[3:6]
        SequenceStructure(sequence='AAA', structure='...')
    """

    sequence: str
    structure: str

    def __post_init__(self) -> None:
        """Validate that sequence and structure have the same length.

        Raises:
            ValueError: If sequence and structure lengths don't match.
        """
        if len(self.sequence) != len(self.structure):
            raise ValueError(
                f"Sequence and structure must have the same length. "
                f"Got sequence length {len(self.sequence)} and structure length {len(self.structure)}."
            )

    def __add__(self, other: "SequenceStructure") -> "SequenceStructure":
        """Concatenate two SequenceStructure objects.

        Args:
            other: Another SequenceStructure to concatenate.

        Returns:
            New SequenceStructure with concatenated sequence and structure.

        Examples:
            >>> ss1 = SequenceStructure("AAA", "...")
            >>> ss2 = SequenceStructure("CCC", "...")
            >>> ss1 + ss2
            SequenceStructure(sequence='AAACCC', structure='......')
        """
        return SequenceStructure(
            self.sequence + other.sequence, self.structure + other.structure
        )

    def __len__(self) -> int:
        """Return the length of the sequence.

        Returns:
            Length of the sequence.
        """
        return len(self.sequence)

    def __getitem__(self, item: int | slice) -> "SequenceStructure":
        """Return a slice of the SequenceStructure.

        Args:
            item: Integer index or slice object.

        Returns:
            New SequenceStructure with sliced sequence and structure.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCC", "(((...)))")
            >>> ss[3:6]
            SequenceStructure(sequence='AAA', structure='...')
        """
        return SequenceStructure(self.sequence[item], self.structure[item])

    def to_dict(self) -> dict[str, str]:
        """Convert to dictionary representation.

        Returns:
            Dictionary with 'sequence' and 'structure' keys.
        """
        return {"sequence": self.sequence, "structure": self.structure}

    def to_comma_delimited(self) -> str:
        """Convert to comma-delimited string representation.

        Returns:
            Comma-separated string: "sequence,structure".

        Examples:
            >>> ss = SequenceStructure("AAA", "...")
            >>> ss.to_comma_delimited()
            'AAA,...'
        """
        return f"{self.sequence},{self.structure}"

    def split_strands(self) -> list["SequenceStructure"]:
        """Split multi-strand structure into individual strands.

        Splits on the strand separator (default: "&") and returns a list
        of SequenceStructure objects, one for each strand.

        Returns:
            List of SequenceStructure objects, one per strand.

        Examples:
            >>> ss = SequenceStructure("GGG&CCC", "(((&)))")
            >>> strands = ss.split_strands()
            >>> len(strands)
            2
        """
        seqs = self.sequence.split(STRAND_SEPARATOR)
        structs = self.structure.split(STRAND_SEPARATOR)
        if len(seqs) != len(structs):
            raise ValueError(
                f"Mismatch between number of sequence strands ({len(seqs)}) "
                f"and structure strands ({len(structs)})."
            )
        return [SequenceStructure(s, st) for s, st in zip(seqs, structs, strict=True)]

    def insert(self, pos: int, other: "SequenceStructure") -> "SequenceStructure":
        """Insert another SequenceStructure at the specified position.

        Args:
            pos: Position to insert at.
            other: SequenceStructure to insert.

        Returns:
            New SequenceStructure with inserted sequence and structure.

        Examples:
            >>> ss1 = SequenceStructure("AAA", "...")
            >>> ss2 = SequenceStructure("GGG", "...")
            >>> ss1.insert(1, ss2)
            SequenceStructure(sequence='AGGGAA', structure='......')
        """
        if pos < 0:
            pos = len(self) + pos
        if pos < 0 or pos > len(self):
            raise IndexError(f"Insert position {pos} out of range [0, {len(self)}]")
        seq = self.sequence[:pos] + other.sequence + self.sequence[pos:]
        struct = self.structure[:pos] + other.structure + self.structure[pos:]
        return SequenceStructure(seq, struct)

    def join(self, other: "SequenceStructure") -> "SequenceStructure":
        """Join two SequenceStructure objects with strand separator.

        Joins the sequences and structures with "&" separator, useful for
        representing multi-strand structures.

        Args:
            other: SequenceStructure to join with.

        Returns:
            New SequenceStructure with joined sequences separated by "&".

        Examples:
            >>> ss1 = SequenceStructure("GGG", "...")
            >>> ss2 = SequenceStructure("CCC", "...")
            >>> ss1.join(ss2)
            SequenceStructure(sequence='GGG&CCC', structure='...&...')
        """
        return SequenceStructure(
            self.sequence + STRAND_SEPARATOR + other.sequence,
            self.structure + STRAND_SEPARATOR + other.structure,
        )

    def replace(self, other: "SequenceStructure", pos: int) -> "SequenceStructure":
        """Replace sequence and structure starting at the specified position.

        Args:
            other: SequenceStructure containing replacement sequence/structure.
            pos: Starting position for replacement.

        Returns:
            New SequenceStructure with replaced sequence and structure.

        Raises:
            ValueError: If position is out of range or replacement extends beyond sequence.

        Examples:
            >>> ss1 = SequenceStructure("ATCG", "....")
            >>> ss2 = SequenceStructure("AA", "()")
            >>> ss1.replace(ss2, 2)
            SequenceStructure(sequence='ATAA', structure='..()')
        """
        if pos < 0:
            pos = len(self) + pos
        if pos < 0 or pos > len(self):
            raise ValueError(
                f"Replacement position {pos} out of range [0, {len(self)}]"
            )
        if pos + len(other) > len(self):
            raise ValueError(
                f"Replacement extends beyond sequence. Position {pos} + length {len(other)} "
                f"> sequence length {len(self)}"
            )
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

    def contains(
        self,
        other: "SequenceStructure",
        start: int | None = None,
        end: int | None = None,
    ) -> bool:
        """Check if this SequenceStructure contains the given substructure.

        Args:
            other: SequenceStructure to search for.
            start: Start position for search (default: 0).
            end: End position for search (default: end of structure).

        Returns:
            True if substructure is found, False otherwise.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCC", "(((...)))")
            >>> sub = SequenceStructure("GAAAC", "(...)")
            >>> ss.contains(sub)
            True
        """
        return len(find(self, other, start, end)) > 0

    def count(
        self,
        other: "SequenceStructure",
        start: int | None = None,
        end: int | None = None,
    ) -> int:
        """Count occurrences of a substructure within this SequenceStructure.

        Args:
            other: SequenceStructure to count.
            start: Start position for search (default: 0).
            end: End position for search (default: end of structure).

        Returns:
            Number of occurrences found.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCCGGGAAACCC", "(((...)))(((...)))")
            >>> sub = SequenceStructure("GAAAC", "(...)")
            >>> ss.count(sub)
            2
        """
        return len(find(self, other, start, end))

    def index(
        self,
        other: "SequenceStructure",
        start: int | None = None,
        end: int | None = None,
    ) -> int:
        """Find the first occurrence of a substructure.

        Args:
            other: SequenceStructure to search for.
            start: Start position for search (default: 0).
            end: End position for search (default: end of structure).

        Returns:
            Starting position of first match.

        Raises:
            ValueError: If substructure is not found.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCC", "(((...)))")
            >>> sub = SequenceStructure("GAAAC", "(...)")
            >>> ss.index(sub)
            2
        """
        matches = find(self, other, start, end)
        if not matches:
            raise ValueError(f"Substructure not found: {other}")
        # For single-strand search, return start position of first match
        if len(matches[0]) == 1:
            return matches[0][0][0]
        # For multi-strand, return start of first strand
        return matches[0][0][0]

    def find_all(
        self,
        other: "SequenceStructure",
        start: int | None = None,
        end: int | None = None,
    ) -> list[tuple[list[int], ...]]:
        """Find all occurrences of a substructure (alias for find function).

        Args:
            other: SequenceStructure to search for.
            start: Start position for search (default: 0).
            end: End position for search (default: end of structure).

        Returns:
            List of matches, same format as find().

        Examples:
            >>> ss = SequenceStructure("GGGAAACCCGGGAAACCC", "(((...)))(((...)))")
            >>> sub = SequenceStructure("GAAAC", "(...)")
            >>> ss.find_all(sub)
            [([2, 7],), ([11, 16],)]
        """
        return find(self, other, start, end)

    def search_by_sequence(
        self, pattern: str, start: int | None = None, end: int | None = None
    ) -> list[tuple[int, int]]:
        """Search for a sequence pattern (ignoring structure).

        Args:
            pattern: Sequence pattern to search for (supports 'N' as wildcard).
            start: Start position for search (default: 0).
            end: End position for search (default: end of structure).

        Returns:
            List of (start, end) positions for each match.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCC", "(((...)))")
            >>> ss.search_by_sequence("AAA")
            [(3, 6)]
            >>> ss.search_by_sequence("NAN")
            [(3, 6)]
        """
        if start is None:
            start = 0
        if end is None:
            end = len(self)
        if start < 0 or end > len(self) or start > end:
            raise ValueError(
                f"Invalid search range: start={start}, end={end}, "
                f"structure length={len(self)}"
            )
        # Handle wildcard N before escaping, then escape special regex characters
        # Use a unique placeholder that won't appear in sequences
        placeholder = "__WILDCARD_N__"
        wildcard_pattern = pattern.replace("N", placeholder)
        escaped_pattern = re.escape(wildcard_pattern).replace(placeholder, r"\S")
        pattern_re = re.compile(escaped_pattern)
        matches = []
        search_seq = self.sequence[start:end]
        for match in pattern_re.finditer(search_seq):
            matches.append((match.start() + start, match.end() + start))
        return matches

    def search_by_structure(
        self, pattern: str, start: int | None = None, end: int | None = None
    ) -> list[tuple[int, int]]:
        """Search for a structure pattern (ignoring sequence).

        Args:
            pattern: Structure pattern to search for (e.g., "(((...)))").
            start: Start position for search (default: 0).
            end: End position for search (default: end of structure).

        Returns:
            List of (start, end) positions for each match.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCC", "(((...)))")
            >>> ss.search_by_structure("(...)")
            [(3, 6)]
        """
        if start is None:
            start = 0
        if end is None:
            end = len(self)
        if start < 0 or end > len(self) or start > end:
            raise ValueError(
                f"Invalid search range: start={start}, end={end}, "
                f"structure length={len(self)}"
            )
        # Escape special regex characters in structure
        escaped_pattern = (
            pattern.replace("(", r"\(").replace(")", r"\)").replace(".", r"\.")
        )
        pattern_re = re.compile(escaped_pattern)
        matches = []
        search_struct = self.structure[start:end]
        for match in pattern_re.finditer(search_struct):
            matches.append((match.start() + start, match.end() + start))
        return matches

    def remove(self, start: int, end: int) -> "SequenceStructure":
        """Remove a subsequence from the SequenceStructure.

        Args:
            start: Start position of region to remove.
            end: End position of region to remove.

        Returns:
            New SequenceStructure with the specified region removed.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCC", "(((...)))")
            >>> ss.remove(3, 6)
            SequenceStructure(sequence='GGGCCC', structure='((()))')
        """
        if start < 0:
            start = len(self) + start
        if end < 0:
            end = len(self) + end
        if start < 0 or end > len(self) or start > end:
            raise ValueError(
                f"Invalid removal range: start={start}, end={end}, "
                f"structure length={len(self)}"
            )
        sequence = self.sequence[:start] + self.sequence[end:]
        structure = self.structure[:start] + self.structure[end:]
        return SequenceStructure(sequence, structure)

    def reverse(self) -> "SequenceStructure":
        """Reverse both sequence and structure.

        Returns:
            New SequenceStructure with reversed sequence and structure.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCC", "(((...)))")
            >>> ss.reverse()
            SequenceStructure(sequence='CCCAAAGGG', structure=')))...(((')
        """
        return SequenceStructure(self.sequence[::-1], self.structure[::-1])

    def strip(self, unpaired_only: bool = True) -> "SequenceStructure":
        """Remove leading and trailing unpaired bases (.) from structure.

        Args:
            unpaired_only: If True, only strip unpaired bases. If False, strip all bases.

        Returns:
            New SequenceStructure with stripped ends.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCC", "..(((...)))..")
            >>> ss.strip()
            SequenceStructure(sequence='GGAAACCC', structure='(((...)))')
        """
        if not unpaired_only:
            # Strip all bases
            return SequenceStructure("", "")
        # Find start: first non-dot character
        start = 0
        while start < len(self) and self.structure[start] == ".":
            start += 1
        # Find end: last non-dot character
        end = len(self)
        while end > start and self.structure[end - 1] == ".":
            end -= 1
        return self[start:end]

    def get_paired_positions(self) -> list[tuple[int, int]]:
        """Get all base pair positions (i, j) where i < j.

        Returns:
            List of (i, j) tuples representing base pairs.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCC", "(((...)))")
            >>> pairs = ss.get_paired_positions()
            >>> len(pairs)
            3
        """
        pairs = []
        stack = []
        for i, char in enumerate(self.structure):
            if char == "(":
                stack.append(i)
            elif char == ")":
                if stack:
                    j = stack.pop()
                    pairs.append((j, i))
        return pairs

    def count_paired_bases(self) -> int:
        """Count the number of paired bases in the structure.

        Returns:
            Number of bases that are part of base pairs.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCC", "(((...)))")
            >>> ss.count_paired_bases()
            6
        """
        return sum(1 for char in self.structure if char in "()")

    def extract_loops(self) -> list["SequenceStructure"]:
        """Extract all hairpin loops from the structure.

        Returns:
            List of SequenceStructure objects, one for each loop.

        Examples:
            >>> ss = SequenceStructure("GGGAAACCC", "(((...)))")
            >>> loops = ss.extract_loops()
            >>> len(loops)
            1
            >>> loops[0].sequence
            'AAA'
        """
        loops = []
        stack = []
        for i, char in enumerate(self.structure):
            if char == "(":
                stack.append(i)
            elif char == ")":
                if stack:
                    start = stack.pop()
                    # Extract the loop region between the opening and closing parentheses
                    loop_start = start + 1
                    loop_end = i
                    # Check if this is a hairpin loop (all unpaired bases)
                    if loop_start < loop_end:
                        loop_seq = self.sequence[loop_start:loop_end]
                        loop_struct = self.structure[loop_start:loop_end]
                        if all(c == "." for c in loop_struct):
                            loops.append(SequenceStructure(loop_seq, loop_struct))
        return loops

    def reverse_complement(self, ntype: str = "DNA") -> "SequenceStructure":
        """Reverse complement the sequence while maintaining structure relationships.

        Args:
            ntype: Type of nucleic acid: "DNA" or "RNA" (default: "DNA").

        Returns:
            New SequenceStructure with reverse complemented sequence and reversed structure.

        Examples:
            >>> ss = SequenceStructure("ATCG", "....")
            >>> ss.reverse_complement()
            SequenceStructure(sequence='CGAT', structure='....')
        """
        # Import here to avoid circular import
        from seq_tools.sequence import get_reverse_complement

        rev_comp_seq = get_reverse_complement(self.sequence, ntype)
        # Reverse the structure to maintain pairing relationships
        rev_struct = self.structure[::-1]
        return SequenceStructure(rev_comp_seq, rev_struct)


def find_seq_struct(
    struct: SequenceStructure,
    sub: SequenceStructure,
    start: int | None = None,
    end: int | None = None,
    format: str = "legacy",
) -> list[tuple[list[int], ...]] | list[Match]:
    """Find positions of a substructure within a structure.

    Alias for `find` function for backward compatibility.

    Args:
        struct: The structure to search within.
        sub: The substructure to search for.
        start: Start position for search (default: 0).
        end: End position for search (default: end of structure).
        format: Output format. "legacy" (default) returns the old format:
                list of tuples with lists. "match" returns a list of Match objects
                which are more readable.

    Returns:
        If format="legacy": List of matches, where each match is a tuple of strand
        positions. Each strand position is [start, end] pair.

        If format="match": List of Match objects with a cleaner interface.

    Examples:
        >>> struct = SequenceStructure("GGGAAACCC", "(((...)))")
        >>> sub = SequenceStructure("GAAAC", "(...)")
        >>> find_seq_struct(struct, sub, format="match")
        [Match(start=2, end=7)]

    See Also:
        find: Main function that performs the search.
    """
    return find(struct, sub, start, end, format=format)


def find(
    struct: SequenceStructure,
    sub: SequenceStructure,
    start: int | None = None,
    end: int | None = None,
    format: str = "legacy",
) -> list[tuple[list[int], ...]] | list[Match]:
    """Find all positions where a substructure appears within a structure.

    Searches for matches where both the sequence pattern and structure pattern
    match. Supports multi-strand searches using "&" separator. Uses regex
    pattern matching with "N" as wildcard in sequences.

    Args:
        struct: The structure to search within.
        sub: The substructure to search for (can be multi-strand).
        start: Start position for search range. If None, starts at 0.
        end: End position for search range. If None, ends at structure length.
        format: Output format. "legacy" (default) returns the old format:
                list of tuples with lists. "match" returns a list of Match objects
                which are more readable.

    Returns:
        If format="legacy": List of all matches. Each match is a tuple of strands,
        where each strand is represented as [start_pos, end_pos]. For single-strand
        searches, returns list of single-element tuples: [([start, end],), ...].
        For multi-strand searches, returns: [([s1_start, s1_end], [s2_start, s2_end], ...), ...].

        If format="match": List of Match objects, which have a cleaner interface.
        Single-strand matches have .start and .end properties, and all matches
        have .strands attribute.

    Examples:
        >>> struct = SequenceStructure("GGGAAACCC", "(((...)))")
        >>> sub = SequenceStructure("GAAAC", "(...)")
        >>> find(struct, sub)
        [[[2, 7]]]
        >>> find(struct, sub, format="match")
        [Match(start=2, end=7)]

        >>> struct = SequenceStructure("GGGAAACCC", "(((...)))")
        >>> sub = SequenceStructure("GGG&CCC", "(((&)))")
        >>> find(struct, sub)
        [[[0, 3], [6, 9]]]
        >>> find(struct, sub, format="match")
        [Match(strands=((0, 3), (6, 9)))]
    """
    if start is None:
        start = 0
    if end is None:
        end = len(struct)

    if start < 0 or end > len(struct) or start > end:
        raise ValueError(
            f"Invalid search range: start={start}, end={end}, "
            f"structure length={len(struct)}"
        )

    if format not in ("legacy", "match"):
        raise ValueError(f'Invalid format: {format}. Must be "legacy" or "match".')

    struct = struct[start:end]
    strands = sub.split_strands()
    strand_matches: list[list[list[int]]] = []

    for strand in strands:
        # Escape special regex characters in sequence, but allow N as wildcard
        placeholder = "__WILDCARD_N__"
        wildcard_seq = strand.sequence.replace("N", placeholder)
        escaped_seq = re.escape(wildcard_seq).replace(placeholder, r"\S")
        pattern_seq = re.compile(r"(?=(" + escaped_seq + r"))")
        # Escape special regex characters in structure
        pattern_ss = re.compile(
            r"(?=("
            + strand.structure.replace("(", r"\(")
            .replace(")", r"\)")
            .replace(".", r"\.")
            + r"))"
        )

        # Find all sequence matches
        matches_seq = [
            str(m.start() + start) + "-" + str(m.start() + len(m.group(1)) + start)
            for m in pattern_seq.finditer(struct.sequence)
        ]

        # Find all structure matches
        matches_ss = [
            str(m.start() + start) + "-" + str(m.start() + len(m.group(1)) + start)
            for m in pattern_ss.finditer(struct.structure)
        ]

        # Find intersection (matches where both sequence and structure match)
        matches = list(set(matches_seq).intersection(set(matches_ss)))
        # Convert string ranges to [start, end] integer pairs
        matches = [list(map(int, i.split("-"))) for i in matches]
        strand_matches.append(matches)

    # Generate all combinations of strand matches using Cartesian product
    # Note: itertools.product returns tuples, preserving original return format
    all_matches = list(itertools.product(*strand_matches))

    if format == "match":
        # Convert to Match objects for cleaner interface
        return [
            Match(strands=tuple((s[0], s[1]) for s in match)) for match in all_matches
        ]

    return all_matches
