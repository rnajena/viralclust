#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
This script mainly exists, because the output of cd-hit-est sucks :)
"""

import sys

import utils

inputFile = sys.argv[1]
sequenceFile = sys.argv[2]

originalHeader = sorted([ header for header, _ in utils.parse_fasta(sequenceFile) ])
truncatedHeader = []
with open(inputFile, 'r') as inputStream:
  for line in inputStream:
    if line.startswith('>'):
      continue
    truncatedHeader.append(line.strip().split('>')[1].split(' ')[0])

truncatedHeader = sorted(truncatedHeader)
assert(len(truncatedHeader) == len(originalHeader))
headerMapping = {x : y for x,y in zip(truncatedHeader,originalHeader)}
assert(len(headerMapping) == len(originalHeader))

with open(inputFile, 'r') as inputSteam:
  for line in inputSteam:
    if line.startswith('>'):
      print(line.strip())
      continue
    accessionID = line.strip().split('>')[1].split(' ')[0]
    line = line.replace(accessionID, headerMapping[accessionID])
    print(line.strip())