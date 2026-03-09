#!/bin/bash

# BSD header for src/ directory
BSD_HEADER="# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: BSD-3-Clause
"

# GPL header for apps/ directory
GPL_HEADER="# SPDX-FileCopyrightText: 2025 Rahil Miten Doshi, Friedrich-Alexander-Universität Erlangen-Nürnberg
# SPDX-License-Identifier: GPL-3.0-or-later
#
# This application depends on waLBerla and pystencils (GPLv3), requiring GPL licensing.
"

# Function to add header to a file
add_header() {
    local file=$1
    local header=$2
    
    # Check if file already has SPDX header
    if head -n 1 "$file" | grep -q "SPDX-FileCopyrightText"; then
        echo "Skipping $file (already has header)"
        return
    fi
    
    # Create temporary file with header + original content
    echo "$header" | cat - "$file" > temp && mv temp "$file"
    echo "Added header to $file"
}

# Export function so it can be used by find -exec
export -f add_header
export BSD_HEADER
export GPL_HEADER

# Add BSD headers to all .py files in src/ (excluding __pycache__)
echo "Adding BSD headers to src/ directory..."
find src/materforge -name "*.py" -type f ! -path "*/__pycache__/*" -exec bash -c 'add_header "$0" "$BSD_HEADER"' {} \;

# Add BSD headers to examples/ directory
echo "Adding BSD headers to examples/ directory..."
find examples -name "*.py" -type f ! -path "*/__pycache__/*" -exec bash -c 'add_header "$0" "$BSD_HEADER"' {} \; 2>/dev/null || true

# Add GPL headers to apps/ directory (excluding walberla submodule)
echo "Adding GPL headers to apps/ directory..."
find apps -name "*.py" -type f ! -path "*/walberla/*" ! -path "*/__pycache__/*" ! -path "*/cmake-build*/*" -exec bash -c 'add_header "$0" "$GPL_HEADER"' {} \;
find apps -name "*.cpp" -type f ! -path "*/walberla/*" ! -path "*/cmake-build*/*" -exec bash -c 'add_header "$0" "$GPL_HEADER"' {} \;

echo "Done! Headers added to all Python files."

