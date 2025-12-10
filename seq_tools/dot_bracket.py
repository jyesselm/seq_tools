"""Simple parsing of dot-bracket notation for RNA secondary structures.

This module provides utilities for converting RNA secondary structure
representations in dot-bracket notation to pair tables.
"""

import collections as col

BRACKET_LEFT = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ"
BRACKET_RIGHT = ")]}>abcdefghijklmnopqrstuvwxyz"


def inverse_brackets(bracket):
    """Create a mapping from bracket characters to their indices.

    Args:
        bracket: String of bracket characters.

    Returns:
        Dictionary mapping each character to its index position.
    """
    res = col.defaultdict(int)
    for i, a in enumerate(bracket):
        res[a] = i
    return res


def dotbracket_to_pairtable(struct):
    """Convert dot-bracket notation to a pair table.

    Converts RNA secondary structure in dot-bracket format to ViennaRNA
    pair table format. Supports multiple bracket types and multi-strand
    structures (with & separator).

    Args:
        struct: Secondary structure in dot-bracket notation.

    Returns:
        Pair table as a list where pt[i] is the index of the paired base,
        or -1 if unpaired.

    Raises:
        ValueError: If structure is empty, has unbalanced brackets, or is invalid.
    """
    if len(struct) == 0:
        raise ValueError("Cannot convert empty structure to pairtable")
    pt = [0] * ((len(struct)) - struct.count("&"))
    # pt[0] = len(struct) - struct.count("&")

    stack = col.defaultdict(list)
    inverse_bracket_left = inverse_brackets(BRACKET_LEFT)
    inverse_bracket_right = inverse_brackets(BRACKET_RIGHT)

    i = 0
    for a in struct:
        if a == "&":
            continue
        i += 1
        # print i,a, pt
        if a == ".":
            pt[i - 1] = -1
        else:
            if a in inverse_bracket_left:
                stack[inverse_bracket_left[a]].append(i)
            else:
                assert a in inverse_bracket_right
                if len(stack[inverse_bracket_right[a]]) == 0:
                    raise ValueError("Too many closing brackets!")
                j = stack[inverse_bracket_right[a]].pop()
                pt[i - 1] = j - 1
                pt[j - 1] = i - 1

    if len(stack[inverse_bracket_left[a]]) != 0:
        raise ValueError("Too many opening brackets!")

    return pt
