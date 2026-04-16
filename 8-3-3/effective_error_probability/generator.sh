#!/bin/bash
cd "$(dirname "$0")"

# Array of custom double values in scientific notation (as strings)
values=("1e-4" "3e-4" "5e-4" "1e-3" "2e-3")

# Loop through the values
for value in "${values[@]}"
do
    # Use awk to convert scientific notation to decimal
    reenc=$(awk -v val="$value" 'BEGIN {printf "%.10f", val}')

    echo "Reencoding error = $reenc"
    
    ./main.exe "1erasure" $reenc
    ./main.exe "flag" $reenc

done