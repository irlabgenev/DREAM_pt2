#!/usr/bin/env python3

'''
USAGE
    fasta_rename.py [OPTION] TABLE [FILE...]

DESCRIPTION
    Replace header in the input FASTA file using a correspondance table.

    TABLE file must list correspondances separated by a tabulation in 
    the following order:
    
    <old name>TAB<new name>
    ...
    
OPTIONS
    --help
        Display this message

'''

import getopt, sys, fileinput, fasta

class Options(dict):

    def __init__(self, argv):
        
        # set default
        self.set_default()
        
        # handle options with getopt
        try:
            opts, args = getopt.getopt(argv[1:], "", ['help'])
        except getopt.GetoptError as e:
            sys.stderr.write(str(e) + '\n' + __doc__)
            sys.exit(1)

        for o, a in opts:
            if o == '--help':
                sys.stdout.write(__doc__)
                sys.exit(0)

        self.args = args
    
    def set_default(self):
    
        # default parameter value
        pass
    
def get_correspondance_table(f):
    return dict( line.strip().split("\t") for line in f )

def main(argv=sys.argv):
    
    # read options and remove options strings from argv (avoid option 
    # names and arguments to be handled as file names by
    # fileinput.input().
    options = Options(argv)
    
    # get conversion table
    with open(options.args.pop(0)) as f:
        d = get_correspondance_table(f)
    sys.argv[1:] = options.args
    
    # organize the main job...
    for header, sequence in fasta.reader(fileinput.input()):
        try:
            new_header = d[header]
        except KeyError:
            new_header = header
        fasta.writer(new_header, sequence, sys.stdout)
    
    # return 0 if everything succeeded
    return 0

# does not execute main if the script is imported as a module
if __name__ == '__main__': sys.exit(main())

