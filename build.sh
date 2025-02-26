#!/bin/bash
set -e

# Install Python and pip
export PYTHON_VERSION=3.9
curl -fsSL https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz | tar xz
cd Python-$PYTHON_VERSION
./configure && make && sudo make install
cd ..

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt

# Build MkDocs site
mkdocs build
