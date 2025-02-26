#!/bin/bash
set -e  # Exit if any command fails

# Ensure Python and pip are installed
if ! command -v python3 &>/dev/null; then
    echo "Python3 is not installed. Installing..."
    apt update && apt install -y python3 python3-pip
fi

# Upgrade pip and install dependencies
python3 -m pip install --upgrade pip
python3 -m pip install mkdocs mkdocs-material mkdocs-awesome-pages-plugin

# Build MkDocs site
mkdocs build

# Ensure output directory exists
mkdir -p site
