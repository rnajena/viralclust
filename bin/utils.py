#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de


import sys
import re


genbankACCRegex = re.compile(r'[^A-Z]([A-Z][0-9]{5}|[A-Z]{2}[0-9]{6}|[A-Z]{2}[0-9]{8}|NC_[0-9]{6})[^0-9]')

def reverseComplement(sequence):
  """

  """
  comp = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A'}
  return(''.join([comp.get(x,'N')  for x in sequence.upper()[::-1]]))



def parse_fasta(filePath):
  """
  """

  with open(filePath, 'r') as inputStream:
    header = ''
    seq = ''

    for line in inputStream:
      if line.startswith(">"):
        if header:
          seq = seq #+ "X"*10 + rev_comp(seq)
          yield (header, seq)

        header = line.rstrip("\n").replace(':','_').replace(' ','_').lstrip(">")
        accessionID = set(re.findall(genbankACCRegex, header))
        if len(accessionID) == 1:
          header = accessionID.pop()
        seq = ''
      else:
        seq += line.rstrip("\n").upper().replace('U','T')

    seq = seq #+ "X"*10 + rev_comp(seq)
    yield (header, seq)