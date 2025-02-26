#!/bin/bash
set -e  # Exit on error

# Install Python manually since Vercel doesn't have it
curl -sS https://bootstrap.pypa.io/get-pip.py | python3

# Ensure Python is installed
if ! command -v python3 &>/dev/null; then
    echo "Python3 is not installed. Installing..."
    apt update && apt install -y python3 python3-pip
fi

# Upgrade pip and install dependencies
python3 -m pip install --upgrade pip
python3 -m pip install mkdocs mkdocs-material mkdocs-awesome-pages-plugin

# Ensure MkDocs is recognized
export PATH=$HOME/.local/bin:$PATH

# Run MkDocs build
python3 -m mkdocs build

# Ensure output directory exists
mkdir -p site
