#!/bin/bash
set -e  # Stop script if any command fails

# Check if Python is installed
if ! command -v python3 &>/dev/null; then
    echo "Python3 is not installed. Installing now..."
    apt update && apt install -y python3 python3-pip
fi

# Upgrade pip & install dependencies
pip3 install --upgrade pip
pip3 install -r requirements.txt

# Build the MkDocs site
mkdocs build

# Ensure output directory exists
mkdir -p site
