#!/usr/bin/env python3

'''
USAGE
    fasta_2nuc.py [OPTION] [FILE...]

DESCRIPTION
    Format a FASTA file into the PAML format (phylip interleaved with 
    names of max length 30).
    
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

def nuc_writer(header, sequence, fout):
    wrap = fasta.textwrap.wrap(sequence, 30, break_on_hyphens=False)
    fout.write("{}\n{}\n".format(header, "\n".join(wrap)))

def main(argv=sys.argv):
    
    # read options and remove options strings from argv (avoid option 
    # names and arguments to be handled as file names by
    # fileinput.input().
    options = Options(argv)
    
    # organize the main job...
    sequences = list(fasta.reader(fileinput.input()))
    n = len(sequences)
    l = len(sequences[0][1])
    
    # check that all sequences are of the same length
    if not all( len(sequence) for header, sequence in sequences ):
        raise ValueError("All sequence are not of the same length.")
    
    # write nuc file header
    sys.stdout.write("   {} {}\n".format(n, l))
    
    # write sequence lines
    for header, sequence in sequences:
        if len(header) > 30:
            raise ValueError("Sequence header length must not exceed 30")
        nuc_writer(header, sequence, sys.stdout)
    
    # return 0 if everything succeeded
    return 0

# does not execute main if the script is imported as a module
if __name__ == '__main__': sys.exit(main())

