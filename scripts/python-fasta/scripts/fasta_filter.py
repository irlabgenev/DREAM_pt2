#!/usr/bin/env python3

'''
USAGE
    fasta_filter.py [OPTION] LIST [FILE...]

DESCRIPTION
    Filter a FASTA files with a given list of sequence headers.
    
OPTIONS
    --short
        Read only the header's leading string before the first space in
        the input FASTA file
        
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
            opts, args = getopt.getopt(argv[1:], "", ['short', 'help'])
        except getopt.GetoptError as e:
            sys.stderr.write(str(e) + '\n' + __doc__)
            sys.exit(1)

        for o, a in opts:
            if o == '--help':
                sys.stdout.write(__doc__)
                sys.exit(0)
            elif o == '--short':
                self['short'] = True

        self.args = args
    
    def set_default(self):
    
        # default parameter value
        self['short'] = False
    
def get_list(f):
    return { line.strip() for line in f }

def main(argv=sys.argv):
    
    # read options and remove options strings from argv (avoid option 
    # names and arguments to be handled as file names by
    # fileinput.input().
    options = Options(argv)
    
    # get conversion table
    with open(options.args.pop(0)) as f:
        s = get_list(f)
    
    sys.argv[1:] = options.args
    
    # organize the main job...
    for header, sequence in fasta.reader(fileinput.input()):
        x = header.split()[0] if options['short'] else header
        if x not in s:
            continue
        fasta.writer(header, sequence, sys.stdout)
    
    # return 0 if everything succeeded
    return 0

# does not execute main if the script is imported as a module
if __name__ == '__main__': sys.exit(main())

