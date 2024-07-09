#!/bin/bash

# Snakemake populates its information as a set of associative arrays.
# These are as follows:
# - snakemake_input
# - snakemake_output
# - snakemake_params
# - snakemake_wildcards
# - snakemake_log
# - snakemake_resources
# - snakemake_config
#
# Note: if you use labels for inputs/outputs/params/etc.,
# there will be two entries for each labeled thing:
# one that's indexed by number, and one that's indexed
# by the label.
#
# Note: nested arrays are not supported, so don't worry about it.

# For convenience, let's capture the first output file as its own variable.
out="${snakemake_output[0]}"

# This is how you get the length of an associative array.
# Notice that each labeled input is counted twice.
echo "Snakemake Input Length:" >> $out
echo ${#snakemake_input[@]} >> $out
echo "" >> $out

# This is how you read all of the keys for an associative array.
# Specifically, it's the "${!snakemake_input[@]}" part, where the "@"
# indicates that Bash should iterate over everything.
# There is always a 0-based index number for each input.
# Labeled inputs are also indexed by their label.
echo "Snakemake Input Keys:" >> $out
for x in "${!snakemake_input[@]}"; do
    echo "$x" >> $out
done
echo "" >> $out

# This is how you read a specific value (by key) from an associative array.
# Specifically, it's the "${snakemake_input[$x]}" part.
# Here, we print each key/value pair.
# Notice again that each labeled input appears twice.
echo "Snakemake Input Values:" >> $out
for x in "${!snakemake_input[@]}"; do
    echo "[$x]=${snakemake_input[$x]}" >> $out
done
