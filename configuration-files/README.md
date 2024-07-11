# Configuration Files
This directory contains examples that show how to use configuration files.
Both YAML and CSV files are used this way.

## Usage
To run this example yourself, simply activate the Conda environment and run Snakemake from this directory.

```shell
# Make sure you've installed the environment,
# and that you're in the "configuration-files" folder.
conda activate snakemake-scratchpad
snakemake --cores 1
```

## Description
Sometimes, you may want to use a separate file to help you specify what a workflow should do.
For example, you may want a separate file that lists the inputs and information about them.
Snakemake supports a few different approaches to this, but we will look at two:
YAML configuration files and Pandas DataFrames (read from a CSV file or similar).

A YAML file is a way to specify key-value information.
This information could be a string, a list, or a set of key-value pairs.
In this way, it works very similar to Python dictionaries and lists.
If you're familiar with JSON, it's very similar to that as well.
In this workflow, we're using it to specify a set of files that need to be merged, the name of the file to produce, and the shell command to use for merging.

Another common approach is to have a CSV file (or something similar) and read it into a Pandas DataFrame.
Snakemake has some useful features for extracting information from a DataFrame.
In this workflow, we're using it to specify lines of text that should go into the files that get merged.
