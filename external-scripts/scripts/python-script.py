# Note: when Snakemake is run, there will be a `snakemake`
# object available without needing to declare it.
# If you're using VSCode (and probably any other IDE), it'll give you
# a warning that `snakemake` hasn't been defined, but that's OK.
#
# The main attributes you may want out of the `snakemake` object are:
# - snakemake.input
# - snakemake.output
# - snakemake.params
# - snakemake.wildcards

def main():
    '''Program execution starts here.'''
    # Open the output file for writing.
    with open(snakemake.output[0], 'w') as f:
        # Snakemake populates all of its information in a `snakemake` object.
        print_attrs(f, snakemake, 'Snakemake Object')
        print_exact(f, type(snakemake), 'Snakemake Object Type')

        # Inputs (along with everything else) are attributes/fields of the object.
        # If you want to ignore any labels, you can simply iterate through
        # `snakemake.input` just like a list.
        print_contents(f, snakemake.input, 'Input List')
        print_attrs(f, snakemake.input, 'Input Object')

        # If you have labels (such as with this example), they are attributes/
        # fields of the `snakemake.input` object.
        print_exact(f, snakemake.input.label1, 'Input Label 1')
        print_contents(f, snakemake.input.label2, 'Input Label 2')
        print_contents(f, snakemake.input.label3, 'Input Label 3')

        # Params (and everything else) work exactly the same as inputs.
        print_contents(f, snakemake.params, 'Params List')
        print_attrs(f, snakemake.params, 'Params Object')
        print_exact(f, snakemake.params.p1, 'Param P1')
        print_exact(f, snakemake.params.p2, 'Param P2')
        print_contents(f, snakemake.params.p3, 'Param P3', skip_newline=True)

def print_contents(f, list, description, skip_newline=False):
    '''
    Helper method for printing the contents of a list,
    skipping any element that starts with an underscore.
    '''
    f.writelines([description, ':\n'])
    f.writelines((f'- {i}\n' for i in list if not str(i).startswith('_')))
    if not skip_newline:
        f.write('\n')

def print_attrs(f, obj, description):
    '''Helper method for printing an object's variables/fields.'''
    print_contents(f, vars(obj), description)

def print_exact(f, obj, description):
    '''Helper method for printing an object itself.'''
    print_contents(f, [obj], description)

# Having a `main` function and calling it this way is not only good
# practice, but also lets you define helper functions after the function
# where it's used, which is nicer to read (in my opinion).
if __name__ == '__main__':
    main()
