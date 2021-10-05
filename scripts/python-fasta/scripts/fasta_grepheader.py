#!/usr/bin/env python3

'''
USAGE
    fasta_grepheader.py [OPTION] PATTERN [FILE...]

DESCRIPTION
    Find sequences using grep pattern to search matching headers.

OPTIONS
    -v, --invert-match
        Select non matching headers

    -x, --line-regexp
        Match full header line
    
    --help
        Display this message

'''

import getopt, sys, fileinput, re, fasta

class Options(dict):

    def __init__(self, argv):
        
        # set default
        self.set_default()
        
        # handle options with getopt
        try:
            opts, args = getopt.getopt(argv[1:], "vx",
            ['invert-match','line-regexp', 'help'])
        except getopt.GetoptError as e:
            sys.stderr.write(str(e) + '\n' + __doc__)
            sys.exit(1)

        for o, a in opts:
            if o == '--help':
                sys.stdout.write(__doc__)
                sys.exit(0)
            elif o in ('-v', '--invert-match'):
                self['inverse'] = True
            elif o in ('-x', '--line-regexp'):
                self['line_regexp'] = True

        self.args = args
    
    def set_default(self):
    
        # default parameter value
        self['inverse'] = False
        self['line_regexp'] = False

def main(argv=sys.argv):
    
    # read options and remove options strings from argv (avoid option 
    # names and arguments to be handled as file names by
    # fileinput.input().
    options = Options(argv)
    if options['line_regexp']:
        p = re.compile('^' + options.args.pop(0) + '$')
    else:
        p = re.compile(options.args.pop(0))
    sys.argv[1:] = options.args
    
    # organize the main job...
    
    # check function
    if options['inverse']:
        check = lambda p, header: p.search(header) is None
    else:
        check = lambda p, header: p.search(header) is not None
    
    # parse input
    for header, sequence in fasta.reader(fileinput.input()):
        if check(p, header):
            fasta.writer(header, sequence, sys.stdout)
    
    # return 0 if everything succeeded
    return 0

# does not execute main if the script is imported as a module
if __name__ == '__main__': sys.exit(main())
