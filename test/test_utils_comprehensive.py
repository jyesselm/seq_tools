"""
Additional comprehensive tests for utils module to improve coverage
"""

from pathlib import Path


# Test the fallback path when importlib.resources is not available
def test_get_resources_path_fallback():
    """Test get_resources_path fallback when files() is not available."""
    import seq_tools.utils as utils_module

    # Save original value
    original_has_files = utils_module._HAS_FILES

    # Temporarily disable files
    utils_module._HAS_FILES = False

    # Import the function after modification
    from seq_tools.utils import get_resources_path

    # Test fallback path
    path = get_resources_path()
    assert isinstance(path, Path)
    assert path.exists()
    assert path.is_dir()

    # Restore original value
    utils_module._HAS_FILES = original_has_files


def test_get_resources_path_with_files():
    """Test get_resources_path with files() available."""
    import seq_tools.utils as utils_module

    # Save original value
    original_has_files = utils_module._HAS_FILES

    # Ensure files is enabled
    utils_module._HAS_FILES = True

    from seq_tools.utils import get_resources_path

    # Test with files
    path = get_resources_path()
    assert isinstance(path, Path)
    assert path.exists()

    # Restore original value
    utils_module._HAS_FILES = original_has_files
