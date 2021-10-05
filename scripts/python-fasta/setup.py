#!/usr/bin/env python3

from distutils.core import setup
import os
setup(name='fasta',
      version='1.0',
      py_modules=['fasta'],
      scripts=[ os.path.join("scripts", x) for x in os.listdir("scripts") ]
      )
