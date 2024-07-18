# Snakemake Scratchpad
A scratchpad for different ways of using [Snakemake](https://snakemake.github.io/), which is a useful tool for data processing workflows (or pipelines) that contain multiple steps.

## Examples in this Repo
Below is a summary of the examples contained in this repo.

- [Basics](./basics/README.md)
  - Some of the basics for using Snakemake.
- [Config Files](./config-files/README.md)
  - Configuring how a workflow runs based on an external file.
  - Both YAML and CSV files are demonstrated.
- [External Scripts](./external-scripts/README.md)
  - Using Python or Bash scripts that are outside the workflow file.
  - R scripts could be used in the same way, but I don't show that here.
  - These scripts show how to retrieve useful information (e.g., inputs and outputs) from the script.
- [Directory Outputs](./directory-outputs/README.md)
  - Using directories as outputs of rule.
  - Referencing directory outputs as inputs to another rule.
  - Also shows how input functions and `checkpoint` rules work.

## Coming Soon
Here is a list of examples that I'll be generating sometime in the near future.

- Using profiles
- Using the SFTP storage plugin
- Using the Slurm executor plugin
- Importing rules from another file
- Generating reports

## Conda Environment
A Conda environment is provided from when this repository was set up.
To create and activate the environment, run these commands (respectively):

```shell
conda env create -f environment.yml
conda activate snakemake-scratchpad
```

If the packages are out of date, run the commands below to update them.

```shell
conda activate snakemake-scratchpad
conda env update --prune
```

If the environment file doesn't work anymore (or is out of date), you can run the following to manually create the environment (and enter "yes" when prompted):

```shell
conda create --name snakemake-scratchpad --channel conda-forge --channel bioconda snakemake pandas
```

If you haven't used Conda before, see the link below for information on how to install it.
I like to use Miniconda, as it installs the minimum necessary to run it via the command line.
If you want to use a GUI, then I'd recommend installing Anaconda.

If you need to install conda, see the link below.
I like to use the Miniconda installation, as it only installs the minimum necessary to run 

https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html
