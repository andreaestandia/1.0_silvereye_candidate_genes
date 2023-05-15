#!/bin/bash

# Usage: ./remove_lines_with_strings.sh strings_file input_file output_file

strings_file="$1"
input_file="$2"
output_file="$3"

# Create a pattern of the strings to be removed
pattern=$(sed 's/.*/\\b&\\b/' "$strings_file" | tr '\n' '|' | sed 's/|$//')

# Remove lines that contain any of the strings
grep -vE "$pattern" "$input_file" > "$output_file"
