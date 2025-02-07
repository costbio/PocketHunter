#!/bin/bash

# Set the installation directory within the repo
INSTALL_DIR="$(pwd)/tools/p2rank"

# Create the installation directory if it doesn't exist
mkdir -p "$INSTALL_DIR"

# Download and extract P2Rank into the repo folder
echo "Downloading P2Rank..."
wget -q -O "$INSTALL_DIR/p2rank_2.5.tar.gz" https://github.com/rdk/p2rank/releases/download/2.5/p2rank_2.5.tar.gz

echo "Extracting P2Rank..."
tar -xzf "$INSTALL_DIR/p2rank_2.5.tar.gz" -C "$INSTALL_DIR" --strip-components=1

# Clean up the tarball
rm "$INSTALL_DIR/p2rank_2.5.tar.gz"

# Ensure prank is executable
chmod +x "$INSTALL_DIR/prank"

# Install Python dependencies if any
if [ -f "requirements.txt" ]; then
    echo "Installing Python dependencies..."
    pip install -r requirements.txt
else
    echo "No Python dependencies to install."
fi

