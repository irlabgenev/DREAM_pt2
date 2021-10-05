#!/usr/bin/env python3

'''
USAGE
    aa_dist.py [OPTION] METHOD [FILE...]

DESCRIPTION
    Compute a pairwise distance matrix between provided sequences using
    one of the different methods. By default, the method used is 
    Levenshtein. For this, sequences do not need to be aligned. Other
    methods may require alignment
    
    Available methods:
    
    Grantham
    Miyata
    Levenshtein
    
OPTIONS
    -d, --delimiter=STR
        Output field delimiter.
    
    -g, --gap
        In amino-acid models that do not incorporate gap, score gap 
        according to the mean substitution score with other amino
        acids. If this option is not selected, gaps are ignored.
    
    --unique-pair-label
        Instead of having the two first columns with the pair labels,
        make a single column with the pair with the names ordered in
        alphanumeric order, separated by a column character.
    
    --help
        Display this message
    
'''

import getopt, sys, fileinput, fasta
from Levenshtein import distance as Levenshtein_distance
from io import StringIO
from functools import partial


MIYATA = StringIO('''aa	Cys	Pro	Ala	Gly	Ser	Thr	Gln	Glu	Asn	Asp	His	Lys	Arg	Val	Leu	Ile	Met	Phe	Tyr	Trp
Cys	0	1.33	1.39	2.22	2.84	1.45	2.48	3.26	2.83	3.48	2.56	3.27	3.06	0.86	1.65	1.63	1.46	2.24	2.38	3.34
Pro	1.33	0	0.06	0.97	0.56	0.87	1.92	2.48	1.8	2.4	2.15	2.94	2.9	1.79	2.7	2.62	2.36	3.17	3.12	4.17
Ala	1.39	0.06	0	0.91	0.51	0.9	1.92	2.46	1.78	2.37	2.17	2.96	2.92	1.85	2.76	2.69	2.42	3.23	3.18	4.23
Gly	2.22	0.97	0.91	0	0.85	1.7	2.48	2.78	1.96	2.37	2.78	3.54	3.58	2.76	3.67	3.6	3.34	4.14	4.08	5.13
Ser	2.84	0.56	0.51	0.85	0	0.89	1.65	2.06	1.31	1.87	1.94	2.71	2.74	2.15	3.04	2.95	2.67	3.45	3.33	4.38
Thr	1.45	0.87	0.9	1.7	0.89	0	1.12	1.83	1.4	2.05	1.32	2.1	2.03	1.42	2.25	2.14	1.86	2.6	2.45	3.5
Gln	2.48	1.92	1.92	2.48	1.65	1.12	0	0.84	0.99	1.47	0.32	1.06	1.13	2.13	2.7	2.57	2.3	2.81	2.48	3.42
Glu	3.26	2.48	2.46	2.78	2.06	1.83	0.84	0	0.85	0.9	0.96	1.14	1.45	2.97	3.53	3.39	3.13	3.59	3.22	4.08
Asn	2.83	1.8	1.78	1.96	1.31	1.4	0.99	0.85	0	0.65	1.29	1.84	2.04	2.76	3.49	3.37	3.08	3.7	3.42	4.39
Asp	3.48	2.4	2.37	2.37	1.87	2.05	1.47	0.9	0.65	0	1.72	2.05	2.34	3.4	4.1	3.98	3.69	4.27	3.95	4.88
His	2.56	2.15	2.17	2.78	1.94	1.32	0.32	0.96	1.29	1.72	0	0.79	0.82	2.11	2.59	2.45	2.19	2.63	2.27	3.16
Lys	3.27	2.94	2.96	3.54	2.71	2.1	1.06	1.14	1.84	2.05	0.79	0	0.4	2.7	2.98	2.84	2.63	2.85	2.42	3.11
Arg	3.06	2.9	2.92	3.58	2.74	2.03	1.13	1.45	2.04	2.34	0.82	0.4	0	2.43	2.62	2.49	2.29	2.47	2.02	2.72
Val	0.86	1.79	1.85	2.76	2.15	1.42	2.13	2.97	2.76	3.4	2.11	2.7	2.43	0	0.91	0.85	0.62	1.43	1.52	2.51
Leu	1.65	2.7	2.76	3.67	3.04	2.25	2.7	3.53	3.49	4.1	2.59	2.98	2.62	0.91	0	0.14	0.41	0.63	0.94	1.73
Ile	1.63	2.62	2.69	3.6	2.95	2.14	2.57	3.39	3.37	3.98	2.45	2.84	2.49	0.85	0.14	0	0.29	0.61	0.86	1.72
Met	1.46	2.36	2.42	3.34	2.67	1.86	2.3	3.13	3.08	3.69	2.19	2.63	2.29	0.62	0.41	0.29	0	0.82	0.93	1.89
Phe	2.24	3.17	3.23	4.14	3.45	2.6	2.81	3.59	3.7	4.27	2.63	2.85	2.47	1.43	0.63	0.61	0.82	0	0.48	1.11
Tyr	2.38	3.12	3.18	4.08	3.33	2.45	2.48	3.22	3.42	3.95	2.27	2.42	2.02	1.52	0.94	0.86	0.93	0.48	0	1.06
Trp	3.34	4.17	4.23	5.13	4.38	3.5	3.42	4.08	4.39	4.88	3.16	3.11	2.72	2.51	1.73	1.72	1.89	1.11	1.06	0
''')

GRANTHAM = StringIO('''aa	Ser	Arg	Leu	Pro	Thr	Ala	Val	Gly	Ile	Phe	Tyr	Cys	His	Gln	Asn	Lys	Asp	Glu	Met	Trp
Ser	0	110	145	74	58	99	124	56	142	155	144	112	89	68	46	121	65	80	135	177
Arg	110	0	102	103	71	112	96	125	97	97	77	180	29	43	86	26	96	54	91	101
Leu	145	102	0	98	92	96	32	138	5	22	36	198	99	113	153	107	172	138	15	61
Pro	74	103	98	0	38	27	68	42	95	114	110	169	77	76	91	103	108	93	87	147
Thr	58	71	92	38	0	58	69	59	89	103	92	149	47	42	65	78	85	65	81	128
Ala	99	112	96	27	58	0	64	60	94	113	112	195	86	91	111	106	126	107	84	148
Val	124	96	32	68	69	64	0	109	29	50	55	192	84	96	133	97	152	121	21	88
Gly	56	125	138	42	59	60	109	0	135	153	147	159	98	87	80	127	94	98	127	184
Ile	142	97	5	95	89	94	29	135	0	21	33	198	94	109	149	102	168	134	10	61
Phe	155	97	22	114	103	113	50	153	21	0	22	205	100	116	158	102	177	140	28	40
Tyr	144	77	36	110	92	112	55	147	33	22	0	194	83	99	143	85	160	122	36	37
Cys	112	180	198	169	149	195	192	159	198	205	194	0	174	154	139	202	154	170	196	215
His	89	29	99	77	47	86	84	98	94	100	83	174	0	24	68	32	81	40	87	115
Gln	68	43	113	76	42	91	96	87	109	116	99	154	24	0	46	53	61	29	101	130
Asn	46	86	153	91	65	111	133	80	149	158	143	139	68	46	0	94	23	42	142	174
Lys	121	26	107	103	78	106	97	127	102	102	85	202	32	53	94	0	101	56	95	110
Asp	65	96	172	108	85	126	152	94	168	177	160	154	81	61	23	101	0	45	160	181
Glu	80	54	138	93	65	107	121	98	134	140	122	170	40	29	42	56	45	0	126	152
Met	135	91	15	87	81	84	21	127	10	28	36	196	87	101	142	95	160	126	0	67
Trp	177	101	61	147	128	148	88	184	61	40	37	215	115	130	174	110	181	152	67	0
''')

aa_names = dict(
         Ala = "A",
         Arg = "R",
         Asn = "N",
         Asp = "D",
         Cys = "C",
         Gln = "Q",
         Glu = "E",
         Gly = "G",
         His = "H",
         Ile = "I",
         Leu = "L",
         Lys = "K",
         Met = "M",
         Phe = "F",
         Pro = "P",
         Ser = "S",
         Thr = "T",
         Trp = "W",
         Tyr = "Y",
         Val = "V")

def make_matrix(M, aa_names=aa_names):
    header = M.readline().split()[1:]
    d = dict( (aa_names[x], dict()) for x in header )
    for line in M:
        cells = line.split()
        aa1 = aa_names[cells[0]]
        scores = map(float, cells[1:])
        for aa2, x in zip(header, scores):
            aa2 = aa_names[aa2]
            d[aa1][aa2] = x
    return d
    
class Options(dict):

    def __init__(self, argv):
        
        # set default
        self.set_default()
        
        # handle options with getopt
        try:
            opts, args = getopt.getopt(argv[1:], "d:g",
                ['delimiter=', 'gap', 'unique-pair-label', 'help'])
        except getopt.GetoptError as e:
            sys.stderr.write(str(e) + '\n' + __doc__)
            sys.exit(1)

        for o, a in opts:
            if o == '--help':
                sys.stdout.write(__doc__)
                sys.exit(0)
            elif o in ('-d', '--delimiter'):
                self["sep"] = a
            elif o in ('-g', '--gap'):
                self["ignore_gaps"] = False
            elif o == "--unique-pair-label":
                self["unique_pair_label"] = True

        self.args = args
    
    def set_default(self):
    
        # default parameter value
        self["ignore_gaps"] = True
        self["sep"] = "\t"
        self["unique_pair_label"] = False

def distance(matrix, seq1, seq2, ignore_gaps=True):
    dist = 0
    
    if not ignore_gaps:
        aa_names = matrix.keys()
        n = len(aa_names)-1
        aa_mean_scores = dict()
        for aa_1 in aa_names:
            s = sum( matrix[aa_1][aa_2] for aa_2 in aa_names if aa_2 != aa_1 )/n
            aa_mean_scores[aa_1] = s
    
    for i in range(len(seq1)):
        if seq1[i] == "-":
            if ignore_gaps or seq2[i] == "-": continue
            dist += aa_mean_scores[seq2[i]]
        elif seq2[i] == "-":
            if ignore_gaps: continue
            dist += aa_mean_scores[seq1[i]]
        else:    
            dist += matrix[seq1[i]][seq2[i]]
    return dist

def main(argv=sys.argv):
    
    # read options and remove options strings from argv (avoid option 
    # names and arguments to be handled as file names by
    # fileinput.input().
    options = Options(argv)
    method = options.args[0]
    sys.argv[1:] = options.args[1:]
    
    global GRANTHAM, MIYATA, distance
    
    # prepare distance function
    if method == "Levenshtein":
        distance = Levenshtein_distance
    elif method == "Grantham":
        distance = partial(distance, 
                           make_matrix(GRANTHAM), 
                           ignore_gaps=options["ignore_gaps"])
    elif method == "Miyata":
        distance = partial(distance, 
                           make_matrix(MIYATA),
                           ignore_gaps=options["ignore_gaps"])
    else:
        sys.stderr.write(f'Unkown method: "{method}"\n')
        return 1
    
    # read sequences
    sequences = list(fasta.reader(fileinput.input()))
    n = len(sequences)
    
    # compute pairwise distances (output long format)
    for i in range(n):
        for j in range(i,n):
            dist = distance(sequences[i][1], sequences[j][1])
            if options["unique_pair_label"]:
                label = ":".join(sorted([sequences[i][0], sequences[j][0]]))
            else:
                label = options['sep'].join([sequences[i][0], sequences[j][0]])
            sys.stdout.write(f"{label}{options['sep']}{dist}\n")
    
    # return 0 if everything succeeded
    return 0

# does not execute main if the script is imported as a module
if __name__ == '__main__': sys.exit(main())

