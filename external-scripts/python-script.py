# Note: when Snakemake is run, there will be a `snakemake`
# object available without needing to declare it.
# The main attributes you may want out of it are:
# - snakemake.input
# - snakemake.output
# - snakemake.params
# - snakemake.wildcards
#
# In this script, we'll just print a bunch of info about
# the inputs and params. Everything else works about the same.

def print_contents(f, list, description):
    f.writelines([description, ':\n'])
    f.writelines((f'- {i}\n' for i in list if not str(i).startswith('_')))
    f.write('\n')

def print_attrs(f, obj, description):
    print_contents(f, dir(obj), description)

def print_exact(f, obj, description):
    print_contents(f, [obj], description)

with open(snakemake.output[0], 'w') as f:
    print_attrs(f, snakemake, 'Snakemake Object')

    print_contents(f, snakemake.input, 'Input List')
    print_attrs(f, snakemake.input, 'Input Object')

    print_exact(f, snakemake.input.label1, 'Input Label 1')
    print_contents(f, snakemake.input.label2, 'Input Label 2')
    print_contents(f, snakemake.input.label3, 'Input Label 3')

    print_contents(f, snakemake.params, 'Params List')
    print_attrs(f, snakemake.params, 'Params Object')

    print_exact(f, snakemake.params.p1, 'Param P1')
    print_exact(f, snakemake.params.p2, 'Param P2')
    print_contents(f, snakemake.params.p3, 'Param P3')
