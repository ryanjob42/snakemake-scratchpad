# Directory Outputs
This directory contains examples that show how to specify entire directories as outputs of a rule in a Snakemake workflow.
These examples also show how to use these directories as inputs to other rules.
It also demonstrates how input functions and `checkpoint` rules work.

## Usage
To run this example yourself, simply activate the Conda environment and run Snakemake from this directory.

```shell
# Make sure you've installed the environment,
# and that you're in the "directory-outputs" folder.
conda activate snakemake-scratchpad
snakemake --cores 1
```

## Description
This workflow has a rule which creates a directory and a random number of files inside that directory (and a subdirectory).
Since we don't know what the outputs will be beforehand, we can only specify the directory as the output to this rule.

After that is run, we have one rule that lists all the contents of the folder and writes it to a file.
This rule takes in the full directory as an input.

At the same time, another rule will merge the contents of all the files into a new file.
This rule uses an input function and relies on the very first rule being a `checkpoint` rule to capture each file as an input.
