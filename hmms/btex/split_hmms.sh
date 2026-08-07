#!/usr/bin/env bash

# Split a concatenated HMMER file into individual .hmm files.
# Each output filename is taken from the model's NAME field.
#
# Usage:
#   bash split_hmms.sh all_models.hmm

set -euo pipefail

input_file="${1:-all_models.hmm}"

if [[ ! -f "$input_file" ]]; then
    echo "Error: file not found: $input_file" >&2
    exit 1
fi

awk '
/^NAME[[:space:]]+/ {
    name = $2
    filename = name ".hmm"
}

{
    if (filename != "") {
        print $0 > filename
    }
}

/^\/\/[[:space:]]*$/ {
    if (filename != "") {
        close(filename)
    }
    filename = ""
    name = ""
}
' "$input_file"

echo "Finished splitting HMMs from: $input_file"
