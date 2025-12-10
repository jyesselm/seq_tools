"""Logger configuration and utilities for seq_tools.

This module provides functions to set up and retrieve loggers with
consistent formatting across the seq_tools package.
"""

import logging
import sys

# logging #####################################################################

APP_LOGGER_NAME = "seq-tools"


def setup_applevel_logger(logger_name=APP_LOGGER_NAME, is_debug=False, file_name=None):
    """Set up an application-level logger with consistent formatting.

    Args:
        logger_name: Name for the logger. Defaults to "seq-tools".
        is_debug: If True, set log level to DEBUG; otherwise INFO. Defaults to False.
        file_name: Optional file path for log output. If None, logs to console only.

    Returns:
        Configured logger instance.
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if is_debug else logging.INFO)

    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    # pylint: disable=C0103
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(sh)

    if file_name:
        # pylint: disable=C0103
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def get_logger(module_name=""):
    """Get or create a logger for a specific module.

    Args:
        module_name: Name of the module requesting the logger. Defaults to empty string.

    Returns:
        Logger instance for the specified module.
    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)
