#!/bin/bash

# Delete build artifacts
rm -rf /path/to/delta-debugger/cppython/build

# Built Cython
cd /path/to/delta-debugger/cppython
python setup.py build_ext --inplace

# Install pip module
cd /path/to/project_root
pip install -e .


