#!/bin/bash
set -e  # Exit on error

# Ensure Python is installed
if ! command -v python3 &>/dev/null; then
    echo "Python3 is not installed. Installing..."
    apt update && apt install -y python3 python3-pip
fi

# Upgrade pip and install dependencies
python3 -m pip install --upgrade pip
python3 -m pip install mkdocs mkdocs-material mkdocs-awesome-pages-plugin

# Find Python's binary path
PYTHON_BIN=$(python3 -c "import sys; print(sys.executable)")

# Build MkDocs site using explicit Python path
$PYTHON_BIN -m mkdocs build

# Ensure output directory exists
mkdir -p site
