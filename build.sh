#!/bin/bash
set -e  # Stop on error

# Install Python & pip (if not installed)
if ! command -v python3 &>/dev/null; then
    echo "Python3 is not installed. Installing..."
    apt update && apt install -y python3 python3-pip
fi

# Upgrade pip and install MkDocs
python3 -m pip install --upgrade pip
python3 -m pip install mkdocs mkdocs-material mkdocs-awesome-pages-plugin

# Ensure MkDocs is accessible
export PATH=$HOME/.local/bin:$PATH

# Run MkDocs build
python3 -m mkdocs build

# Ensure output directory exists
mkdir -p site
