#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.rst').read()
doclink = """
Documentation
-------------

The full documentation is at http://seq_tools.rtfd.org."""
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='seq_tools',
    version='0.1.0',
    description='simple functions for manipulating sequences and secondary structures in pandas dataframe format',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Joe Yesselman',
    author_email='jyesselm@unl.edu',
    url='https://github.com/jyesselm/seq_tools',
    packages=[
        'seq_tools',
    ],
    package_dir={'seq_tools': 'seq_tools'},
    py_modules=[
        'seq_tools/data_frame', 'seq_tools/dot_bracket', 'seq_tools/extinction_coeff',
        'seq_tools/seq_tools', 'seq_tools/sequence'
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='seq_tools',
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    entry_points = {
        'console_scripts' : [
            'seq_tools = seq_tools.seq_tools:main',
        ]
    }
)
