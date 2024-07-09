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

out="${snakemake_output[0]}"

echo "Snakemake Input Length:" >> $out
echo ${#snakemake_input[@]} >> $out
echo "" >> $out

echo "Snakemake Input Keys:" >> $out
for x in "${!snakemake_input[@]}"; do
    echo "$x" >> $out
done
echo "" >> $out

echo "Snakemake Input Values:" >> $out
for x in "${!snakemake_input[@]}"; do
    echo "[$x]=${snakemake_input[$x]}" >> $out
done
