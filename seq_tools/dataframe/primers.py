"""Primer and common sequence handling for dataframes."""

import pandas as pd
import vienna

from seq_tools.config import T7_PROMOTER
from seq_tools.dataframe.core import fold, run_in_parallel, trim
from seq_tools.dataframe.transforms import to_rna
from seq_tools.utils import get_resources_path
from seq_tools.validation import validate_dataframe


def has_5p_sequence(df: pd.DataFrame, p5_seq: str) -> bool:
    """Check if all sequences start with the given 5' sequence.

    Args:
        df: DataFrame containing sequences.
        p5_seq: 5' sequence to check for.

    Returns:
        True if all sequences start with p5_seq, False otherwise.
    """
    return df["sequence"].str.startswith(p5_seq).all()


def has_3p_sequence(df: pd.DataFrame, p3_seq: str) -> bool:
    """Check if all sequences end with the given 3' sequence.

    Args:
        df: DataFrame containing sequences.
        p3_seq: 3' sequence to check for.

    Returns:
        True if all sequences end with p3_seq, False otherwise.
    """
    return df["sequence"].str.endswith(p3_seq).all()


def has_sequence(df: pd.DataFrame, seq: str) -> bool:
    """Check if all sequences contain the given subsequence.

    Args:
        df: DataFrame containing sequences.
        seq: Subsequence to search for.

    Returns:
        True if all sequences contain seq, False otherwise.
    """
    return df["sequence"].str.contains(seq).all()


def has_t7_promoter(df: pd.DataFrame) -> bool:
    """Check if all sequences in the dataframe have a T7 promoter.

    Args:
        df: DataFrame with sequences.

    Returns:
        True if all sequences start with T7 promoter, False otherwise.
    """
    has_t7 = df[df["sequence"].str.startswith(T7_PROMOTER)]
    if len(has_t7) != len(df):
        return False
    return True


def find_longest_common_prefix(df: pd.DataFrame) -> str:
    """Find the longest common prefix sequence from all sequences in the DataFrame.

    Args:
        df: DataFrame containing sequences.

    Returns:
        The longest common prefix sequence, or empty string if none found.
    """
    if len(df) == 0:
        return ""

    sequences = df["sequence"].tolist()
    if len(sequences) == 0:
        return ""

    if len(sequences) == 1:
        return sequences[0]

    # Find the minimum length
    min_len = min(len(seq) for seq in sequences)

    # Check each position
    common_prefix = ""
    for i in range(min_len):
        # Check if all sequences have the same character at position i
        char = sequences[0][i]
        if all(seq[i] == char for seq in sequences):
            common_prefix += char
        else:
            break

    return common_prefix


def find_longest_common_suffix(df: pd.DataFrame) -> str:
    """Find the longest common suffix sequence from all sequences in the DataFrame.

    Args:
        df: DataFrame containing sequences.

    Returns:
        The longest common suffix sequence, or empty string if none found.
    """
    if len(df) == 0:
        return ""

    sequences = df["sequence"].tolist()
    if len(sequences) == 0:
        return ""

    if len(sequences) == 1:
        return sequences[0]

    # Find the minimum length
    min_len = min(len(seq) for seq in sequences)

    # Check each position from the end
    common_suffix = ""
    for i in range(1, min_len + 1):
        # Check if all sequences have the same character at position -i
        char = sequences[0][-i]
        if all(seq[-i] == char for seq in sequences):
            common_suffix = char + common_suffix
        else:
            break

    return common_suffix


def trim_p5_and_p3(
    df: pd.DataFrame, extra_columns: list[str] | None = None
) -> pd.DataFrame:
    """Trim the 5' and 3' ends of the data in the DataFrame.

    This function reads a CSV file containing p5 sequences, converts these
    sequences to RNA, checks for a common p5 sequence in the given DataFrame,
    and trims the DataFrame based on the length of this common p5 sequence and
    a fixed 3' end length.

    Args:
        df: A DataFrame with a 'data' column containing sequences as strings.
        extra_columns: Additional columns to trim along with sequence.

    Returns:
        A trimmed DataFrame with the 5' and 3' ends trimmed.

    Raises:
        ValueError: If no common p5 sequence is found or the sequence is not
                    registered in the CSV file.
    """
    if extra_columns is None:
        extra_columns = []
    df_p5 = pd.read_csv(get_resources_path() / "p5_sequences.csv")
    df_p3 = pd.read_csv(get_resources_path() / "p3_sequences.csv")
    common_p5_seq = ""
    common_p3_seq = ""
    for p5_seq in df_p5["sequence"]:
        if has_5p_sequence(df, p5_seq):
            common_p5_seq = p5_seq
    for p3_seq in df_p3["sequence"]:
        if has_3p_sequence(df, p3_seq):
            common_p3_seq = p3_seq
    if len(common_p5_seq) == 0 or len(common_p3_seq) == 0:
        raise ValueError("No common p5 or p3 sequence found")
    return trim(df, len(common_p5_seq), len(common_p3_seq), extra_columns)


def remove_common_p5_p3(
    df: pd.DataFrame, extra_columns: list[str] | None = None
) -> pd.DataFrame:
    """Identify and remove common 5' and 3' sequences from the DataFrame.

    This function reads CSV files containing p5 and p3 sequences, identifies
    which common sequences are present in all sequences in the DataFrame,
    and removes them.

    Args:
        df: DataFrame containing sequences.
        extra_columns: Additional columns to trim along with sequence.

    Returns:
        DataFrame with common 5' and 3' sequences removed.

    Raises:
        ValueError: If no common p5 or p3 sequence is found.
    """
    if extra_columns is None:
        extra_columns = []
    df_p5 = pd.read_csv(get_resources_path() / "p5_sequences.csv")
    df_p3 = pd.read_csv(get_resources_path() / "p3_sequences.csv")

    common_p5_seq = ""
    common_p3_seq = ""

    # Find common p5 sequence (check all sequences)
    for _, row in df_p5.iterrows():
        p5_seq = row["sequence"]
        if has_5p_sequence(df, p5_seq):
            common_p5_seq = p5_seq
            break

    # Find common p3 sequence (check all sequences)
    for _, row in df_p3.iterrows():
        p3_seq = row["sequence"]
        if has_3p_sequence(df, p3_seq):
            common_p3_seq = p3_seq
            break

    if len(common_p5_seq) == 0 and len(common_p3_seq) == 0:
        raise ValueError("No common p5 or p3 sequence found")

    return trim(df, len(common_p5_seq), len(common_p3_seq), extra_columns)


def _check_p5_match(df: pd.DataFrame, p5_seq: str, p5_struct: str) -> tuple[bool, bool]:
    """Check if p5 sequence and structure match.

    Args:
        df: DataFrame with sequences and optionally structures.
        p5_seq: 5' sequence to check.
        p5_struct: 5' structure to check (can be empty).

    Returns:
        Tuple of (sequence_match, structure_match).
    """
    has_structure = "structure" in df.columns

    # Check sequence match
    if not has_5p_sequence(df, p5_seq):
        return False, False

    seq_match = True
    struct_match = False

    # Check structure match if both are available
    if p5_struct and has_structure:
        p5_struct_len = len(p5_struct)
        if p5_struct_len > 0:
            struct_match = df["structure"].str[:p5_struct_len] == p5_struct
            struct_match = struct_match.all()

    return seq_match, struct_match


def _check_p3_match(df: pd.DataFrame, p3_seq: str, p3_struct: str) -> tuple[bool, bool]:
    """Check if p3 sequence and structure match.

    Args:
        df: DataFrame with sequences and optionally structures.
        p3_seq: 3' sequence to check.
        p3_struct: 3' structure to check (can be empty).

    Returns:
        Tuple of (sequence_match, structure_match).
    """
    has_structure = "structure" in df.columns

    # Check sequence match
    if not has_3p_sequence(df, p3_seq):
        return False, False

    seq_match = True
    struct_match = False

    # Check structure match if both are available
    if p3_struct and has_structure:
        p3_struct_len = len(p3_struct)
        if p3_struct_len > 0:
            struct_match = df["structure"].str[-p3_struct_len:] == p3_struct
            struct_match = struct_match.all()

    return seq_match, struct_match


def remove_common_p5_p3_by_structure(
    df: pd.DataFrame, extra_columns: list[str] | None = None
) -> pd.DataFrame:
    """Identify and remove common 5' and 3' sequences based on both sequence and structure.

    This function matches sequences that have both matching sequence and structure
    patterns from the resource files. It checks the p5_sequences.csv which contains
    both sequence and structure information.

    Args:
        df: DataFrame containing sequences and optionally structures.
        extra_columns: Additional columns to trim along with sequence.

    Returns:
        DataFrame with common 5' and 3' sequences removed based on structure matching.

    Raises:
        ValueError: If no matching sequence/structure pattern is found, or if
                    structure column is missing when required.
    """
    if extra_columns is None:
        extra_columns = []
    df_p5 = pd.read_csv(get_resources_path() / "p5_sequences.csv")
    df_p3 = pd.read_csv(get_resources_path() / "p3_sequences.csv")

    common_p5_seq = ""
    common_p3_seq = ""

    # Find common p5 sequence and structure
    for _, row in df_p5.iterrows():
        p5_seq = row["sequence"]
        p5_struct = row.get("structure", "")

        seq_match, struct_match = _check_p5_match(df, p5_seq, p5_struct)

        if seq_match and (struct_match or not p5_struct):
            common_p5_seq = p5_seq
            break

    # Find common p3 sequence and structure
    for _, row in df_p3.iterrows():
        p3_seq = row["sequence"]
        p3_struct = row.get("structure", "")

        seq_match, struct_match = _check_p3_match(df, p3_seq, p3_struct)

        if seq_match and (struct_match or not p3_struct):
            common_p3_seq = p3_seq
            break

    if len(common_p5_seq) == 0 and len(common_p3_seq) == 0:
        raise ValueError("No common p5 or p3 sequence/structure pattern found")

    return trim(df, len(common_p5_seq), len(common_p3_seq), extra_columns)


def _generate_match_warnings(
    seq: str,
    struct: str,
    is_p5: bool,
    seq_match: bool,
    struct_match: bool,
    has_structure: bool,
) -> list[str]:
    """Generate warnings for sequence/structure matching.

    Args:
        seq: The sequence being matched.
        struct: The structure being matched.
        is_p5: True if this is a 5' sequence, False for 3'.
        seq_match: Whether sequence matched.
        struct_match: Whether structure matched.
        has_structure: Whether dataframe has structure column.

    Returns:
        List of warning messages.
    """
    warnings = []
    end = "5'" if is_p5 else "3'"

    if seq_match and not struct_match:
        if struct and has_structure:
            warnings.append(
                f"{end} sequence '{seq}' matches but structure doesn't match"
            )
        elif struct and not has_structure:
            warnings.append(
                f"{end} sequence '{seq}' matches but no structure column found in input"
            )

    return warnings


def remove_common_seqs(
    df: pd.DataFrame, extra_columns: list[str] | None = None
) -> tuple[pd.DataFrame, dict]:
    """Identify and remove common 5' and 3' sequences, checking both sequence and structure.

    This function attempts to match both sequence and structure patterns. It provides
    warnings if only sequence or only structure matches, but will still proceed with
    removal based on sequence match if structure is not available or doesn't match.

    Args:
        df: DataFrame containing sequences and optionally structures.
        extra_columns: Additional columns to trim along with sequence.

    Returns:
        Tuple of (DataFrame with common sequences removed, dict with warnings).

    Raises:
        ValueError: If no common p5 or p3 sequence is found.
    """
    if extra_columns is None:
        extra_columns = []
    df_p5 = pd.read_csv(get_resources_path() / "p5_sequences.csv")
    df_p3 = pd.read_csv(get_resources_path() / "p3_sequences.csv")

    common_p5_seq = ""
    common_p3_seq = ""
    warnings = []
    has_structure = "structure" in df.columns

    # Find common p5 sequence and structure
    for _, row in df_p5.iterrows():
        p5_seq = row["sequence"]
        p5_struct = row.get("structure", "")

        seq_match, struct_match = _check_p5_match(df, p5_seq, p5_struct)

        if seq_match:
            common_p5_seq = p5_seq
            warnings.extend(
                _generate_match_warnings(
                    p5_seq, p5_struct, True, seq_match, struct_match, has_structure
                )
            )
            break

    # Find common p3 sequence and structure
    for _, row in df_p3.iterrows():
        p3_seq = row["sequence"]
        p3_struct = row.get("structure", "")

        seq_match, struct_match = _check_p3_match(df, p3_seq, p3_struct)

        if seq_match:
            common_p3_seq = p3_seq
            warnings.extend(
                _generate_match_warnings(
                    p3_seq, p3_struct, False, seq_match, struct_match, has_structure
                )
            )
            break

    if len(common_p5_seq) == 0 and len(common_p3_seq) == 0:
        raise ValueError("No common p5 or p3 sequence found")

    removed_info = {}
    if common_p5_seq:
        removed_info["p5_sequence"] = common_p5_seq
        removed_info["p5_length"] = len(common_p5_seq)
    if common_p3_seq:
        removed_info["p3_sequence"] = common_p3_seq
        removed_info["p3_length"] = len(common_p3_seq)

    return trim(df, len(common_p5_seq), len(common_p3_seq), extra_columns), {
        "warnings": warnings,
        "removed": removed_info,
    }


def transcribe(df: pd.DataFrame, ignore_missing_t7: bool = False) -> pd.DataFrame:
    """Transcribe DNA template sequences to RNA and remove T7 promoter.

    Converts DNA to RNA, removes the T7 promoter (first 20 bases), and folds
    the resulting RNA sequences.

    Args:
        df: DataFrame with DNA template sequences.
        ignore_missing_t7: If True, proceed even if sequences don't have T7 promoter.

    Returns:
        DataFrame with transcribed RNA sequences and folded structures.

    Raises:
        ValueError: If not all sequences have T7 promoter and ignore_missing_t7 is False.
    """
    if not has_t7_promoter(df) and not ignore_missing_t7:
        raise ValueError("not all sequences start with T7 promoter")
    if not ignore_missing_t7:
        df = trim(df, 20, 0)
    df = to_rna(df)
    df = fold(df)
    return df


def _fold_and_validate(row: pd.Series) -> pd.Series:
    """Fold and validate a single row.

    Args:
        row: DataFrame row with sequence and expected_structure.

    Returns:
        Series with predicted_structure, mfe, ens_defect, and structure_match.
    """
    v_res = vienna.fold(row["sequence"])
    predicted_structure = v_res.dot_bracket
    structure_match = predicted_structure == row["expected_structure"]
    return pd.Series(
        {
            "predicted_structure": predicted_structure,
            "mfe": v_res.mfe,
            "ens_defect": v_res.ens_defect,
            "structure_match": structure_match,
        }
    )


def add_common_seqs(
    df: pd.DataFrame,
    p5_seq: str,
    p5_structure: str,
    p3_seq: str,
    p3_structure: str,
    parallel: bool = False,
    workers: int = 4,
) -> pd.DataFrame:
    """Add common 5' and/or 3' sequences to sequences and validate structure.

    This function adds the specified p5 and/or p3 sequences to each sequence in the
    dataframe, folds the resulting RNA sequences, and validates that the predicted
    structure matches the expected structure (p5_structure + structure + p3_structure).

    Args:
        df: DataFrame with 'sequence' and 'structure' columns.
        p5_seq: 5' sequence to add (can be empty string).
        p5_structure: Expected structure for p5 sequence (can be empty string).
        p3_seq: 3' sequence to add (can be empty string).
        p3_structure: Expected structure for p3 sequence (can be empty string).
        parallel: If True, use parallel processing for structure validation.
        workers: Number of workers to use for parallel processing.

    Returns:
        DataFrame with columns: name, sequence, structure, predicted_structure, structure_match.
    """
    validate_dataframe(df)
    if "structure" not in df.columns:
        raise ValueError("DataFrame must contain a 'structure' column")

    df_result = df.copy()

    # Build expected structure for each row
    df_result["expected_structure"] = df_result["structure"].apply(
        lambda s: p5_structure + s + p3_structure
    )

    # Add sequences to create new sequences
    df_result["sequence"] = df_result["sequence"].apply(lambda x: p5_seq + x + p3_seq)

    if parallel:

        def process_chunk(chunk_df: pd.DataFrame) -> pd.DataFrame:
            result = chunk_df.apply(_fold_and_validate, axis=1)
            return pd.concat([chunk_df, result], axis=1)

        df_result = run_in_parallel(df_result, process_chunk, workers)
    else:
        validation_results = df_result.apply(_fold_and_validate, axis=1)
        df_result = pd.concat([df_result, validation_results], axis=1)

    # Clean up: remove expected_structure column as it's internal
    df_result = df_result.drop(columns=["expected_structure"])

    # Reorder columns for better output
    columns_order = [
        "name",
        "sequence",
        "structure",
        "predicted_structure",
        "structure_match",
    ]
    # Add any remaining columns
    for col in df_result.columns:
        if col not in columns_order:
            columns_order.append(col)

    df_result = df_result[columns_order]

    return df_result
