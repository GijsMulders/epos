#!/bin/sh
# creates the documentation using sphinx from the docstrings
sphinx-apidoc -f -o _docs EPOS/
make html
# documentation in _build/html/index.html'