# Coding Standards Quick Reference

## The Rules

### 1. Indentation: Max 3 Levels
```python
# ❌ BAD - 4 levels
def process(data):
    if data:
        for item in data:
            if item.valid:
                if item.score > 10:
                    return True

# ✅ GOOD - 2 levels with early returns
def process(data):
    if not data:
        return False

    for item in data:
        if not item.valid:
            continue
        if item.score > 10:
            return True
    return False
```

### 2. Function Length: ≤ 30 Lines
```python
# ❌ BAD - Too long, does too much
def process_and_validate_and_save(df, p5, p3, validate, parallel, output):
    # ... 50 lines of code ...
    pass

# ✅ GOOD - Split into smaller functions
def add_primers(df: pd.DataFrame, p5: str, p3: str) -> pd.DataFrame:
    """Add primers to sequences."""
    # ... 10 lines ...

def validate_structures(df: pd.DataFrame, parallel: bool) -> pd.DataFrame:
    """Validate RNA structures."""
    # ... 15 lines ...

def save_results(df: pd.DataFrame, output: str) -> None:
    """Save results to file."""
    # ... 5 lines ...
```

### 3. Single Responsibility
```python
# ❌ BAD - Does multiple things
def process_file(filename):
    df = pd.read_csv(filename)
    df = validate(df)
    df = transform(df)
    df.to_csv("output.csv")
    return df

# ✅ GOOD - Each function does one thing
def read_file(filename: str) -> pd.DataFrame:
    """Read sequences from CSV file."""
    return pd.read_csv(filename)

def validate_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """Validate sequence data."""
    # validation logic
    return df

def save_results(df: pd.DataFrame, output: str) -> None:
    """Save results to CSV file."""
    df.to_csv(output, index=False)
```

### 4. Type Hints Required
```python
# ❌ BAD - No type hints
def fold(seq):
    return vienna.fold(seq)

# ✅ GOOD - Clear type hints
def fold_sequence(seq: str) -> tuple[str, float]:
    """Fold RNA sequence and return structure and MFE.

    Args:
        seq: RNA sequence to fold.

    Returns:
        Tuple of (structure, mfe).
    """
    result = vienna.fold(seq)
    return (result.dot_bracket, result.mfe)
```

### 5. Docstrings Required
```python
# ❌ BAD - No docstring
def calc(x, y):
    return x * y + 5

# ✅ GOOD - Clear docstring
def calculate_score(count: int, weight: float) -> float:
    """Calculate weighted score with baseline offset.

    Args:
        count: Number of occurrences.
        weight: Weight factor to apply.

    Returns:
        Weighted score with +5 baseline offset.
    """
    return count * weight + 5
```

### 6. Simple Over Clever
```python
# ❌ BAD - Too clever
result = [x for x in [y.upper() for y in data if y] if len(x) > 3]

# ✅ GOOD - Clear and simple
def filter_valid_sequences(data: list[str]) -> list[str]:
    """Filter sequences that are non-empty and longer than 3."""
    valid = []
    for seq in data:
        if seq and len(seq) > 3:
            valid.append(seq.upper())
    return valid
```

### 7. Module Size: 200-300 Lines Target
```
# ❌ BAD - Giant file
seq_tools/
├── dataframe.py (1200 lines)

# ✅ GOOD - Organized into modules
seq_tools/
├── dataframe/
│   ├── __init__.py
│   ├── folding.py (210 lines)
│   ├── sequences.py (245 lines)
│   ├── validation.py (180 lines)
│   └── generation.py (220 lines)
```

## Pre-Commit Checklist

Run these before every commit:

```bash
make format      # Auto-format with black
make lint        # Check with ruff
make type-check  # Verify types with mypy
make coverage    # Run tests (90% required)
```

Or run all at once:
```bash
make check
```

## Common Patterns

### Early Returns
```python
# ✅ GOOD - Reduces nesting
def validate(df: pd.DataFrame) -> bool:
    if df is None:
        return False
    if "sequence" not in df.columns:
        return False
    if len(df) == 0:
        return False
    return True
```

### Guard Clauses
```python
# ✅ GOOD - Handle edge cases early
def process(data: list[str]) -> list[str]:
    if not data:
        return []

    # Main logic here
    return [d.upper() for d in data]
```

### Extract Functions
```python
# ✅ GOOD - Complex logic in separate function
def validate_structure(predicted: str, expected: str) -> bool:
    """Check if predicted structure matches expected."""
    return predicted == expected

def process_sequences(df: pd.DataFrame) -> pd.DataFrame:
    df["match"] = df.apply(
        lambda row: validate_structure(row["pred"], row["exp"]),
        axis=1
    )
    return df
```

## Tools Configuration

All tools are configured in `pyproject.toml`:
- **Ruff**: Formatting and linting, line length 88, max complexity 10
- **Mypy**: Strict mode, all functions need types
- **Pytest**: 90% coverage minimum

See [CONTRIBUTING.md](../CONTRIBUTING.md) for complete details.
