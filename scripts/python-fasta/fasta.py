#!/usr/bin/env python3

import textwrap, re

SEQ = re.compile("\S+")

def reader(f):
    header, sequence = "", ""
    for line in f:
        if line.startswith(">"): 
            if header and sequence:
                yield header, sequence
                sequence = ""
            header = line.strip()[1:]
        else:
            sequence += "".join(SEQ.findall(line))
    yield header, sequence

def writer(header, sequence, fout):
    wrap = textwrap.wrap(sequence, 60, break_on_hyphens=False)
    fout.write(">{}\n{}\n".format(header, "\n".join(wrap)))

def read_range(s):
    a, b = s.split("-")
    return slice(int(a)-1, int(b))
