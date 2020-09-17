#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

import sys

import utils

sequenceFile = sys.argv[1]
outputFile = sequenceFile.replace("positive","negative")

with open(outputFile, 'w') as outputStream:
  for header, sequence in utils.parse_fasta(sequenceFile):
    outputStream.write(f">{header}\n{utils.reverseComplement(sequence)}\n")
