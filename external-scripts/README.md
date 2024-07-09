# External Scripts
This directory contains examples that show how to use external Bash and Python scripts in a Snakemake workflow.
Other kinds of scripts, such as R or Julia, are also supported by Snakemake, but are not shown here.
For the official documentation of external scripts, see: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts

There are four main rules in the Snakemake file, each of which simply writes information to the output file about the data structures populated by Snakemake.
Three of these main rules are for Bash scripts, and demonstrate some potentially unexpected behavior in how the associative arrays are populated.
The other main rule is for Python.

The `scripts` folder contains the Bash and Python scripts which run the test.
Their outputs are in the `results` folder, and are named similarly to the Snakemake rule that generated them.

The workflow also contains an `all` rule, which runs everything, and a `generate_input` rule, which just creates temporary input files.
