#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# Install Python dependencies
pip install --upgrade pip
pip install -r requirements.txt

# Build MkDocs site
mkdocs build

# Ensure output directory exists
mkdir -p site
