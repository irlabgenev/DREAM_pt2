#!/usr/bin/env python3

'''
USAGE
    fasta_dealign.py [OPTION] [FILE...]

DESCRIPTION
    Remove gaps from FASTA sequences.

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

def main(argv=sys.argv):
    
    # read options and remove options strings from argv (avoid option 
    # names and arguments to be handled as file names by
    # fileinput.input().
    options = Options(argv)
    sys.argv[1:] = options.args
    
    # organize the main job...
    for header, sequence in fasta.reader(fileinput.input()):
        sequence = sequence.replace("-", "")
        fasta.writer(header, sequence, sys.stdout)
    
    # return 0 if everything succeeded
    return 0

# does not execute main if the script is imported as a module
if __name__ == '__main__': sys.exit(main())

