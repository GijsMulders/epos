#!/bin/sh
# This script creates the documentation using autosphinx from the docstrings
# You can access the documentation in _build/html/index.html
sphinx-apidoc -f -o docs EPOS/
make html
