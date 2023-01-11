"""
setup script for seq_tools
"""
# !/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

with open("README.md", "r", encoding="utf-8") as f:
    readme = f.read()

with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = f.read().splitlines()

setup(
    name="rna_seq_tools",
    version="0.6.0",
    description="simple functions for manipulating sequences and secondary "
    "structures in pandas dataframe format",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    url="https://github.com/jyesselm/seq_tools",
    packages=[
        "seq_tools",
    ],
    package_dir={"seq_tools": "seq_tools"},
    py_modules=[
        "seq_tools/dataframe",
        "seq_tools/dot_bracket",
        "seq_tools/cli",
        "seq_tools/extinction_coeff",
        "seq_tools/logger",
        "seq_tools/sequence",
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords="seq_tools",
    classifiers=[
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],
    entry_points={
        "console_scripts": [
            "seq-tools = seq_tools.cli:cli",
        ]
    },
)
