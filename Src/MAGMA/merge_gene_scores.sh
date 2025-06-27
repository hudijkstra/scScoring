#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <directory> <output_file>"
    exit 1
fi

DIR="$1"
OUTPUT_FILE="$2"

# Ensure output file is empty or created
> "$OUTPUT_FILE"

first_file=true

for file in "$DIR"/*genes.out; do
    if [ -f "$file" ]; then
        if $first_file; then
            # Keep header from the first file only
            head -n 1 "$file" >> "$OUTPUT_FILE"
            first_file=false
        fi
        # Skip the header line and append the rest
        tail -n +2 "$file" >> "$OUTPUT_FILE"
    fi
done
