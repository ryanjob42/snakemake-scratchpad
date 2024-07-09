# Snakemake Scratchpad
A scratchpad for different ways of using [Snakemake](https://snakemake.github.io/), which is a useful tool for data processing workflows (or pipelines) that contain multiple steps.

Below is a summary of the tests contained in this repo.

- [External Scripts](./external-scripts/README.md)
  - Demonstrates how to use external scripts.
  - Shows how Snakemake populates useful information (e.g., inputs and outputs) for your script to use.

## Conda Environment
A Conda environment is provided from when this repository was set up.
To create and activate the environment, run these commands (respectively):

```shell
conda env create -f environment.yml
conda activate snakemake-scratchpad
```

If the environment file doesn't work anymore, you can run the following to manually create the environment (and enter "yes" when prompted):

```shell
conda create --name snakemake-scratchpad --channel conda-forge --channel bioconda snakemake
```

If you haven't used Conda before, see the link below for information on how to install it.
I like to use Miniconda, as it installs the minimum necessary to run it via the command line.
If you want to use a GUI, then I'd recommend installing Anaconda.

If you need to install conda, see the link below.
I like to use the Miniconda installation, as it only installs the minimum necessary to run 

https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html
